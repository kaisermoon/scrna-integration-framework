"""PR1c5b: runner action control and best-effort failure audit."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scrna_integration.run_contract import atomic_write_json

ROOT = Path(__file__).resolve().parents[1]


def _runner():
    spec = importlib.util.spec_from_file_location("pr1c5_control", ROOT / "scripts/smoke_run_notebooks.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _setup(monkeypatch, tmp_path: Path):
    runner = _runner()
    (tmp_path / "stage.ipynb").write_text("{}")
    monkeypatch.setattr(runner, "PROJECT_ROOT", str(tmp_path))
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path / "executed"))
    monkeypatch.setattr(runner, "RUN_ROOT", tmp_path / "runs")
    return runner


def _executed(runner, *, error=False, count=1, evalue="boom"):
    output = Path(runner.OUTPUT_DIR) / "stage_exec_stage.ipynb"
    output.parent.mkdir(parents=True, exist_ok=True)
    errors = ([{"output_type": "error", "ename": "RuntimeError", "evalue": evalue,
                "traceback": ["long traceback", "bounded traceback summary"]}] if error else [])
    cell = {"cell_type": "code", "source": ["x = 1"], "execution_count": 1, "outputs": errors}
    output.write_text(json.dumps({"cells": [cell for _ in range(count)]}))


def _partial(root: Path, name="run", *, nested=False):
    draft = root / name / "draft"
    target = draft / "nested" / "partial.bin" if nested else draft / "partial.bin"
    target.parent.mkdir(parents=True)
    target.write_bytes(b"keep")
    return target


def test_partial_cell_error_audit_keeps_details_and_original_error(monkeypatch, tmp_path: Path) -> None:
    runner = _setup(monkeypatch, tmp_path)

    def fake_run(*args, **kwargs):
        _executed(runner, error=True, count=runner.MAX_FAILED_CELLS + 3, evalue="x" * 100_000)
        _partial(runner.RUN_ROOT)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    monkeypatch.setattr(runner, "NOTEBOOKS", [
        {"name": "stage", "notebook": "stage.ipynb", "expected_stage": "04_embedded"}])
    assert runner.main() == 1
    result = json.loads((Path(runner.OUTPUT_DIR) / "smoke_run_summary.json").read_text())["results"][0]
    manifest = json.loads(Path(result["manifest"]).read_text())
    detail = manifest["failure"]["cell_errors"][0]
    assert result["status"] == "ERROR" and runner.TRUNCATION_MARKER in result["error"]
    assert len(detail["evalue"]) <= runner.MAX_EVALUE_CHARS and detail["evalue"].endswith(runner.TRUNCATION_MARKER)
    assert result["failures"] == manifest["failure"]["cell_errors"] and len(result["failures"]) == 20
    assert result["truncated_cell_error_count"] == manifest["failure"]["truncated_cell_error_count"] == 3
    assert list(Path(result["manifest"]).parent.iterdir()) == [Path(result["manifest"])] and not hasattr(runner, "DOWNSTREAM")


@pytest.mark.parametrize("mode", ["write", "cleanup"])
def test_audit_failure_never_replaces_original_summary(monkeypatch, tmp_path: Path, mode) -> None:
    runner = _setup(monkeypatch, tmp_path)
    marker = None

    def fake_run(*args, **kwargs):
        nonlocal marker
        marker = _partial(runner.RUN_ROOT, nested=True)
        return SimpleNamespace(returncode=7, stdout="", stderr="original nbconvert failure")

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    monkeypatch.setattr(runner, "NOTEBOOKS", [
        {"name": "stage", "notebook": "stage.ipynb", "expected_stage": "04_embedded"}])
    target, name, message = ((runner, "atomic_write_json", "disk full") if mode == "write" else
                             (runner.shutil, "rmtree", "busy"))
    monkeypatch.setattr(target, name, lambda *args, **kwargs: (_ for _ in ()).throw(OSError(message)))
    assert runner.main() == 1
    summary = json.loads((Path(runner.OUTPUT_DIR) / "smoke_run_summary.json").read_text())["results"][0]
    assert summary["error"] == "original nbconvert failure" and summary["audit_error"]
    assert marker.read_bytes() == b"keep"
    manifest_path = marker.parents[1] / "manifest.json"
    if mode == "write":
        assert summary["manifest"] is None and not manifest_path.exists()
    else:
        assert summary["manifest"] == str(manifest_path)
        assert json.loads(manifest_path.read_text())["audit_cleanup_errors"]


@pytest.mark.parametrize("mode", ["multiple", "existing", "unsafe", "state-link"])
def test_ambiguous_or_unsafe_partial_draft_is_unchanged(monkeypatch, tmp_path: Path, mode) -> None:
    runner = _setup(monkeypatch, tmp_path)
    markers = []

    def fake_run(*args, **kwargs):
        if mode == "multiple":
            markers.extend([_partial(runner.RUN_ROOT, "one"), _partial(runner.RUN_ROOT, "two")])
        elif mode == "unsafe":
            markers.append(_partial(runner.RUN_ROOT, "bad id"))
        elif mode == "state-link":
            marker = _partial(tmp_path / "outside")
            run_dir = runner.RUN_ROOT / "run"
            run_dir.mkdir(parents=True)
            (run_dir / "draft").symlink_to(marker.parent, target_is_directory=True)
            markers.append(marker)
        else:
            marker = _partial(runner.RUN_ROOT)
            atomic_write_json(marker.parent / "manifest.json",
                              {"run_id": "run", "stage_status": "FAILED", "failure": {"source": "notebook"}})
            markers.append(marker)
        return SimpleNamespace(returncode=4, stdout="", stderr="original")

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    result = runner.run_notebook("stage", "stage.ipynb", "04_embedded")
    assert result["error"] == "original" and all(path.read_bytes() == b"keep" for path in markers)
    if mode == "existing":
        assert json.loads((markers[0].parent / "manifest.json").read_text())["failure"] == {"source": "notebook"}
