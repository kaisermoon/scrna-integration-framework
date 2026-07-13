from __future__ import annotations

import json
from pathlib import Path

import pytest

from scrna_integration.run_contract import (
    MethodStatus,
    StageStatus,
    atomic_write_json,
    determine_stage_status,
    prepare_run,
    promote_run,
    resume_run,
    sha256_file,
    validate_checkpoint,
)

ACCEPTANCE = {"accepted_by": "researcher", "accepted_at": "2026-07-13T12:00:00Z"}


@pytest.mark.parametrize(
    ("required", "hard", "needs_review", "warnings", "expected"),
    [
        ({"pca": MethodStatus.SUCCESS}, {"counts_valid": True}, False, (), StageStatus.SUCCESS),
        ({"pca": MethodStatus.FAILED}, {"counts_valid": True}, False, (), StageStatus.FAILED),
        ({"pca": MethodStatus.UNAVAILABLE}, {"counts_valid": True}, False, (), StageStatus.FAILED),
        ({"pca": MethodStatus.SUCCESS}, {"counts_valid": False}, False, (), StageStatus.FAILED),
        (
            {"pca": MethodStatus.SUCCESS},
            {"counts_valid": True},
            True,
            (),
            StageStatus.NEEDS_REVIEW,
        ),
        (
            {"pca": MethodStatus.SUCCESS},
            {"counts_valid": True},
            False,
            ("optional method failed",),
            StageStatus.SUCCESS_WITH_WARNINGS,
        ),
    ],
)
def test_determine_stage_status(required, hard, needs_review, warnings, expected) -> None:
    actual = determine_stage_status(required, hard, needs_review=needs_review, warnings=warnings)
    assert actual is expected


def test_determine_stage_status_rejects_empty_hard_postconditions() -> None:
    with pytest.raises(ValueError, match="must not be empty"):
        determine_stage_status({}, {})


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


@pytest.mark.parametrize(
    "acceptance",
    [
        {},
        {"accepted_by": "", "accepted_at": "2026-07-13T12:00:00Z"},
        {"accepted_by": "researcher", "accepted_at": ""},
    ],
)
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
        promote_run(
            paths,
            warning_acceptance=ACCEPTANCE,
        )

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
