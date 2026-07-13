"""PR1c5a: core runner accepts only a fresh, validated run manifest."""

from __future__ import annotations

import importlib.util
import inspect
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scrna_integration.run_contract import atomic_write_json, prepare_run, promote_run, sha256_file

ROOT = Path(__file__).resolve().parents[1]


def _runner():
    spec = importlib.util.spec_from_file_location("pr1c5_runner", ROOT / "scripts/smoke_run_notebooks.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _write_run(root: Path, run_id: str, *, stage="05_clustered", status="SUCCESS",
               state="promoted", accepted=True, artifact=False):
    paths = prepare_run(root, run_id)
    checkpoint = paths.draft_dir / "checkpoint.h5ad"
    checkpoint.write_bytes(b"checkpoint")
    manifest = {"run_id": run_id, "stage": stage, "stage_status": status,
                "checkpoint": {"path": checkpoint.name, "sha256": sha256_file(checkpoint)}}
    if status == "SUCCESS_WITH_WARNINGS" and accepted:
        manifest["warning_acceptance"] = {"accepted_by": "PI", "accepted_at": "2026-07-14"}
    if artifact:
        report = paths.draft_dir / "report.json"
        report.write_bytes(b"report")
        manifest["artifacts"] = [{"role": "report", "path": report.name, "sha256": sha256_file(report)}]
    atomic_write_json(paths.manifest_path, manifest)
    if state == "promoted":
        promote_run(paths)
        return paths.promoted_dir / "manifest.json"
    return paths.manifest_path


def _rewrite(path: Path, **changes):
    payload = json.loads(path.read_text())
    payload.update(changes)
    atomic_write_json(path, payload, overwrite=True)


def _error(runner, before, root, code="INVALID_MANIFEST"):
    with pytest.raises(runner.ManifestValidationError) as caught:
        runner.validate_manifest_delta(before, root, "05_clustered")
    assert caught.value.code == code


def test_run_notebook_returns_fresh_valid_promoted_manifest(monkeypatch, tmp_path: Path) -> None:
    runner, notebook = _runner(), tmp_path / "stage.ipynb"
    notebook.write_text("{}")
    monkeypatch.setattr(runner, "PROJECT_ROOT", str(tmp_path))
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path / "executed"))
    monkeypatch.setattr(runner, "RUN_ROOT", tmp_path / "runs")

    def fake_run(*args, **kwargs):
        output = Path(runner.OUTPUT_DIR) / "stage_exec_stage.ipynb"
        output.parent.mkdir(parents=True, exist_ok=True)
        cell = {"cell_type": "code", "source": ["x = 1"], "execution_count": 1, "outputs": []}
        output.write_text(json.dumps({"cells": [cell]}))
        _write_run(runner.RUN_ROOT, "fresh")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    outcome = runner.run_notebook("stage", "stage.ipynb", "05_clustered")["run_outcome"]
    assert [outcome[key] for key in ("run_id", "stage", "state", "status")] == [
        "fresh", "05_clustered", "promoted", "SUCCESS"]
    assert Path(outcome["manifest"]).is_file() and Path(outcome["checkpoint"]).is_file()
    source = inspect.getsource(runner.run_notebook)
    assert all(set(item) == {"name", "notebook", "expected_stage"} for item in runner.NOTEBOOKS)
    assert ".h5ad" not in repr(runner.NOTEBOOKS)
    assert all(token not in source for token in ("st_mtime", "output_exists", "MISSING_OUTPUT", "STALE_OUTPUT"))


def test_legacy_flat_h5ad_and_zero_manifest_are_rejected(tmp_path: Path) -> None:
    runner, root = _runner(), tmp_path / "runs"
    before = runner.snapshot_run_manifests(root)
    (tmp_path / "05_clustered_v1.h5ad").write_bytes(b"legacy")
    _error(runner, before, root, "MISSING_MANIFEST")
    _write_run(root, "one")
    _write_run(root, "two")
    _error(runner, {}, root, "MULTIPLE_MANIFESTS")


def test_old_manifest_unchanged_modified_and_deleted(tmp_path: Path) -> None:
    runner, root = _runner(), tmp_path / "runs"
    old = _write_run(root, "old")
    before = runner.snapshot_run_manifests(root)
    _write_run(root, "fresh")
    assert runner.validate_manifest_delta(before, root, "05_clustered")["run_id"] == "fresh"
    before = runner.snapshot_run_manifests(root)
    old.unlink()
    _error(runner, before, root, "STALE_RUN")
    old = _write_run(root, "another")
    before = runner.snapshot_run_manifests(root)
    _rewrite(old, stage="modified")
    _error(runner, before, root, "STALE_RUN")


@pytest.mark.parametrize("fault", ["stage", "manifest-run-id", "state"])
def test_manifest_identity_and_path_must_match(fault, tmp_path: Path) -> None:
    runner, root = _runner(), tmp_path / "runs"
    if fault == "state":
        path = root / "run" / "review" / "manifest.json"
        path.parent.mkdir(parents=True)
        path.write_text("{}")
    else:
        path = _write_run(root, "run", stage="wrong" if fault == "stage" else "05_clustered")
        if fault == "manifest-run-id":
            _rewrite(path, run_id="other")
    _error(runner, {}, root)


@pytest.mark.parametrize("unsafe", ["bad id", "运行", "bad\nid"])
def test_run_id_directory_requires_safe_ascii(unsafe, tmp_path: Path) -> None:
    runner, root = _runner(), tmp_path / "runs"
    path = _write_run(root, "safe")
    path.parents[1].rename(root / unsafe)
    _error(runner, {}, root)


@pytest.mark.parametrize("status", [None, "UNKNOWN", "FAILED"])
def test_missing_unknown_and_failed_status_are_rejected(status, tmp_path: Path) -> None:
    runner, root = _runner(), tmp_path / "runs"
    path = _write_run(root, "run")
    _rewrite(path, stage_status=status)
    _error(runner, {}, root)


@pytest.mark.parametrize("status", ["SUCCESS", "SUCCESS_WITH_WARNINGS"])
def test_success_states_must_be_promoted(status, tmp_path: Path) -> None:
    runner, root = _runner(), tmp_path / "runs"
    _write_run(root, "draft", status=status, state="draft")
    _error(runner, {}, root)


@pytest.mark.parametrize("target", ["checkpoint", "artifact"])
def test_tampered_checkpoint_or_artifact_is_rejected(target, tmp_path: Path) -> None:
    runner, root = _runner(), tmp_path / "runs"
    path = _write_run(root, "run", artifact=True)
    payload = json.loads(path.read_text())
    relative = payload["checkpoint"]["path"] if target == "checkpoint" else payload["artifacts"][0]["path"]
    (path.parent / relative).write_bytes(b"tampered")
    _error(runner, {}, root)


def test_warning_acceptance_and_needs_review(tmp_path: Path) -> None:
    runner, accepted = _runner(), tmp_path / "accepted"
    _write_run(accepted, "warnings", status="SUCCESS_WITH_WARNINGS")
    assert runner.validate_manifest_delta({}, accepted, "05_clustered")["status"] == "SUCCESS_WITH_WARNINGS"
    rejected = tmp_path / "rejected"
    _write_run(rejected, "warnings", status="SUCCESS_WITH_WARNINGS", state="draft", accepted=False)
    _error(runner, {}, rejected)
    review = tmp_path / "review"
    _write_run(review, "review", status="NEEDS_REVIEW", state="draft")
    assert runner.validate_manifest_delta({}, review, "05_clustered")["state"] == "draft"


@pytest.mark.parametrize("kind", ["state", "manifest"])
def test_state_and_manifest_symlinks_are_rejected(kind, tmp_path: Path) -> None:
    runner, root = _runner(), tmp_path / "runs"
    if kind == "state":
        source = _write_run(tmp_path / "outside", "run")
        (root / "run").mkdir(parents=True)
        (root / "run" / "promoted").symlink_to(source.parent, target_is_directory=True)
    else:
        path = _write_run(root, "run")
        target = tmp_path / "manifest.json"
        path.replace(target)
        path.symlink_to(target)
    _error(runner, {}, root)


def test_manifest_mutation_during_validation_is_stale(monkeypatch, tmp_path: Path) -> None:
    runner, root = _runner(), tmp_path / "runs"
    path = _write_run(root, "run")
    real_validate = runner.validate_checkpoint

    def mutate(manifest_path):
        checkpoint = real_validate(manifest_path)
        _rewrite(path, race=True)
        return checkpoint

    monkeypatch.setattr(runner, "validate_checkpoint", mutate)
    _error(runner, {}, root, "STALE_RUN")
