from __future__ import annotations

import json
import os
import subprocess
import sys
from enum import Enum
from pathlib import Path

import pytest

from scrna_integration.run_contract import (
    MethodStatus,
    StageStatus,
    atomic_write_json,
    collect_runtime_provenance,
    determine_stage_status,
    prepare_run,
    promote_run,
    publish_compatibility_symlink,
    resume_run,
    sha256_file,
    snapshot_effective_parameters,
    validate_checkpoint,
)

ACCEPTANCE = {"accepted_by": "researcher", "accepted_at": "2026-07-13T12:00:00Z"}


class _Choice(Enum):
    FIRST = "first"


class _Scalar:
    def item(self):
        return 3


def test_snapshot_effective_parameters_normalizes_and_sorts(tmp_path: Path) -> None:
    namespace = {
        "ZETA": (_Choice.FIRST, tmp_path / "input.h5ad"),
        "ALPHA": {"nested": [_Scalar(), True]},
        "lowercase": "ignored",
        "RUNTIME": os,
    }

    assert snapshot_effective_parameters(namespace) == {
        "ALPHA": {"nested": [3, True]},
        "ZETA": ["first", str(tmp_path / "input.h5ad")],
    }


def test_snapshot_effective_parameters_redacts_secrets_and_honors_exclude() -> None:
    params = snapshot_effective_parameters(
        {
            "OPENAI_API_KEY": object(),
            "API_TOKEN": "token-value",
            "ACCESS_TOKEN": "token-value",
            "PASSWORD_FILE": Path("secret.txt"),
            "CREDENTIAL": {"raw": "secret"},
            "BATCH_KEY": "batch",
            "LABEL_KEY": "cell_type",
            "GROUP_KEY": "sample",
            "VISIBLE": 1,
            "OMITTED": 2,
        },
        exclude=("OMITTED",),
    )

    assert params == {
        "ACCESS_TOKEN": "<redacted>",
        "API_TOKEN": "<redacted>",
        "BATCH_KEY": "batch",
        "CREDENTIAL": "<redacted>",
        "GROUP_KEY": "sample",
        "LABEL_KEY": "cell_type",
        "OPENAI_API_KEY": "<redacted>",
        "PASSWORD_FILE": "<redacted>",
        "VISIBLE": 1,
    }


def test_snapshot_effective_parameters_rejects_unsupported_value() -> None:
    with pytest.raises(TypeError, match="UNSUPPORTED"):
        snapshot_effective_parameters({"UNSUPPORTED": {1, 2}})


def test_collect_runtime_provenance_records_git_and_packages(tmp_path: Path) -> None:
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "config", "user.email", "test@example.org"], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "config", "user.name", "Test"], check=True)
    (tmp_path / "tracked.txt").write_text("content", encoding="utf-8")
    subprocess.run(["git", "-C", str(tmp_path), "add", "tracked.txt"], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "commit", "-qm", "test"], check=True)
    expected_commit = subprocess.check_output(
        ["git", "-C", str(tmp_path), "rev-parse", "HEAD"], text=True
    ).strip()

    provenance = collect_runtime_provenance(tmp_path, ("pytest", "definitely-missing-package"))

    assert provenance["python"] == sys.version.split()[0]
    assert provenance["platform"]
    assert provenance["git_commit"] == expected_commit
    assert provenance["packages"]["pytest"] != "unavailable"
    assert provenance["packages"]["definitely-missing-package"] == "unavailable"


def test_collect_runtime_provenance_marks_git_unavailable(tmp_path: Path) -> None:
    assert collect_runtime_provenance(tmp_path)["git_commit"] == "unavailable"


def test_publish_compatibility_symlink_creates_and_updates_relative_link(tmp_path: Path) -> None:
    target_one = tmp_path / "runs" / "one.h5ad"
    target_two = tmp_path / "runs" / "two.h5ad"
    target_one.parent.mkdir()
    target_one.write_text("one", encoding="utf-8")
    target_two.write_text("two", encoding="utf-8")
    link = tmp_path / "latest.h5ad"

    publish_compatibility_symlink(link, target_one)
    assert link.is_symlink()
    assert not os.readlink(link).startswith("/")
    assert link.resolve() == target_one.resolve()

    publish_compatibility_symlink(link, target_two)
    assert link.resolve() == target_two.resolve()
    assert list(tmp_path.glob(".latest.h5ad.*.tmp")) == []


def test_publish_compatibility_symlink_rejects_real_file_and_missing_target(
    tmp_path: Path,
) -> None:
    target = tmp_path / "target.h5ad"
    target.write_text("target", encoding="utf-8")
    real_file = tmp_path / "latest.h5ad"
    real_file.write_text("do not replace", encoding="utf-8")

    with pytest.raises(FileExistsError):
        publish_compatibility_symlink(real_file, target)
    with pytest.raises(FileNotFoundError):
        publish_compatibility_symlink(tmp_path / "missing-link", tmp_path / "missing-target")

    assert real_file.read_text(encoding="utf-8") == "do not replace"


def test_publish_compatibility_symlink_cleans_temp_after_replace_failure(
    monkeypatch, tmp_path: Path
) -> None:
    target = tmp_path / "target.h5ad"
    target.write_text("target", encoding="utf-8")

    def fail_replace(source, destination):
        raise OSError("simulated replace failure")

    monkeypatch.setattr(os, "replace", fail_replace)
    with pytest.raises(OSError, match="simulated replace failure"):
        publish_compatibility_symlink(tmp_path / "latest.h5ad", target)

    assert list(tmp_path.glob(".latest.h5ad.*.tmp")) == []


@pytest.mark.parametrize(
    ("required", "hard", "needs_review", "warnings", "expected"),
    [
        ({"pca": MethodStatus.SUCCESS}, {"counts_valid": True}, False, (), StageStatus.SUCCESS),
        ({"pca": MethodStatus.FAILED}, {"counts_valid": True}, False, (), StageStatus.FAILED),
        ({"pca": MethodStatus.UNAVAILABLE}, {"counts_valid": True}, False, (), StageStatus.FAILED),
        ({"pca": MethodStatus.SUCCESS}, {"counts_valid": False}, False, (), StageStatus.FAILED),
        ({"pca": MethodStatus.SUCCESS}, {"counts_valid": True}, True, (), StageStatus.NEEDS_REVIEW),
        ({"pca": MethodStatus.SUCCESS}, {"counts_valid": True}, False, ("warning",), StageStatus.SUCCESS_WITH_WARNINGS),
    ],
)
def test_determine_stage_status(required, hard, needs_review, warnings, expected) -> None:
    actual = determine_stage_status(required, hard, needs_review=needs_review, warnings=warnings)
    assert actual is expected


def test_determine_stage_status_rejects_empty_hard_postconditions() -> None:
    with pytest.raises(ValueError, match="must not be empty"):
        determine_stage_status({}, {})


def test_no_required_methods_must_be_explicitly_allowed() -> None:
    with pytest.raises(ValueError, match="required_methods must not be empty"):
        determine_stage_status({}, {"counts_valid": True})
    assert determine_stage_status(
        {}, {"counts_valid": True}, allow_no_required_methods=True
    ) is StageStatus.SUCCESS


def test_prepare_run_rejects_duplicate_run_id(tmp_path: Path) -> None:
    paths = prepare_run(tmp_path, "run-001")
    assert paths.draft_dir.is_dir()

    with pytest.raises(FileExistsError):
        prepare_run(tmp_path, "run-001")


def test_atomic_json_write_is_no_clobber_and_leaves_no_temp_file(tmp_path: Path) -> None:
    destination = tmp_path / "manifest.json"
    atomic_write_json(destination, {"value": 1})

    with pytest.raises(FileExistsError):
        atomic_write_json(destination, {"value": 2})

    assert json.loads(destination.read_text(encoding="utf-8")) == {"value": 1}
    assert list(tmp_path.glob(".*.tmp")) == []


def _write_draft(tmp_path: Path, status: StageStatus, *, good_hash: bool = True):
    paths = prepare_run(tmp_path, status.value.lower())
    checkpoint = paths.draft_dir / "checkpoint.h5ad"
    checkpoint.write_bytes(b"checkpoint")
    checksum = sha256_file(checkpoint) if good_hash else "0" * 64
    atomic_write_json(
        paths.manifest_path,
        {
            "run_id": paths.run_id,
            "stage_status": status.value,
            "checkpoint": {"path": checkpoint.name, "sha256": checksum},
        },
    )
    return paths, checkpoint


@pytest.mark.parametrize("status", [StageStatus.FAILED, StageStatus.NEEDS_REVIEW])
def test_failed_or_unreviewed_run_cannot_promote(tmp_path: Path, status: StageStatus) -> None:
    paths, _ = _write_draft(tmp_path, status)

    with pytest.raises(ValueError, match="cannot be promoted"):
        promote_run(paths)

    assert not paths.promoted_dir.exists()


def test_warning_run_requires_recorded_acceptance(tmp_path: Path) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS_WITH_WARNINGS)

    with pytest.raises(ValueError, match="requires accepted_by and accepted_at"):
        promote_run(paths)

    promoted = promote_run(paths, warning_acceptance=ACCEPTANCE)
    manifest = json.loads((paths.promoted_dir / "manifest.json").read_text(encoding="utf-8"))
    assert promoted == paths.promoted_dir / "checkpoint.h5ad"
    assert manifest["warning_acceptance"] == ACCEPTANCE


def test_warning_run_reuses_valid_manifest_acceptance(tmp_path: Path) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS_WITH_WARNINGS)
    manifest = json.loads(paths.manifest_path.read_text(encoding="utf-8"))
    manifest["warning_acceptance"] = ACCEPTANCE
    atomic_write_json(paths.manifest_path, manifest, overwrite=True)

    assert promote_run(paths).is_file()


@pytest.mark.parametrize("acceptance", [
    {},
    {"accepted_by": "", "accepted_at": "2026-07-13T12:00:00Z"},
    {"accepted_by": "researcher", "accepted_at": ""},
])
def test_warning_run_rejects_incomplete_acceptance(tmp_path: Path, acceptance) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS_WITH_WARNINGS)

    with pytest.raises(ValueError, match="requires accepted_by and accepted_at"):
        promote_run(paths, warning_acceptance=acceptance)


def test_checkpoint_hash_mismatch_rejects_validation_and_promotion(tmp_path: Path) -> None:
    paths, checkpoint = _write_draft(tmp_path, StageStatus.SUCCESS, good_hash=False)

    with pytest.raises(ValueError, match="hash does not match"):
        validate_checkpoint(paths.manifest_path, checkpoint)
    with pytest.raises(ValueError, match="hash does not match"):
        promote_run(paths)


def test_bad_hash_does_not_record_warning_acceptance(tmp_path: Path) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS_WITH_WARNINGS, good_hash=False)
    before = paths.manifest_path.read_bytes()

    with pytest.raises(ValueError, match="hash does not match"):
        promote_run(paths, warning_acceptance=ACCEPTANCE)

    assert paths.manifest_path.read_bytes() == before


def test_manifest_checkpoint_path_mismatch_is_rejected(tmp_path: Path) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    other = paths.draft_dir / "other.h5ad"
    other.write_bytes(b"checkpoint")

    with pytest.raises(ValueError, match="does not match manifest"):
        validate_checkpoint(paths.manifest_path, other)


def test_resume_run_requires_a_validated_checkpoint(tmp_path: Path) -> None:
    paths, checkpoint = _write_draft(tmp_path, StageStatus.SUCCESS)

    assert resume_run(tmp_path, paths.run_id) == paths

    checkpoint.write_bytes(b"changed")
    with pytest.raises(ValueError, match="hash does not match"):
        resume_run(tmp_path, paths.run_id)


def test_successful_promotion_moves_whole_draft_atomically(tmp_path: Path) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)

    promoted = promote_run(paths)

    assert promoted.read_bytes() == b"checkpoint"
    assert not paths.draft_dir.exists()
    assert paths.promoted_dir.is_dir()
    assert resume_run(tmp_path, paths.run_id, promoted=True) == paths


def test_relative_root_prepare_promote_and_resume(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.chdir(tmp_path)
    paths, _ = _write_draft(Path("results/runs"), StageStatus.SUCCESS)

    promoted = promote_run(paths)

    assert promoted == tmp_path / "results/runs/success/promoted/checkpoint.h5ad"
    assert resume_run("results/runs", paths.run_id, promoted=True) == paths
