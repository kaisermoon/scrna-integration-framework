from __future__ import annotations

import hashlib
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
    aggregate_method_status,
    atomic_write_json,
    collect_runtime_provenance,
    determine_stage_status,
    prepare_run,
    promote_run,
    publish_compatibility_symlink,
    resume_run,
    sha256_file,
    snapshot_effective_parameters,
    validate_artifacts,
    validate_checkpoint,
    validate_expression_contract,
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

    assert snapshot_effective_parameters(namespace, path_root=tmp_path) == {
        "ALPHA": {"nested": [3, True]},
        "ZETA": ["first", "input.h5ad"],
    }


def test_snapshot_effective_parameters_redacts_secrets_and_honors_exclude() -> None:
    params = snapshot_effective_parameters(
        {
            "OPENAI_API_KEY": object(),
            "OPENAI_KEY": "secret",
            "API_TOKEN": "token-value",
            "PRIVATE_KEY": "secret",
            "AWS_ACCESS_KEY_ID": "secret",
            "AUTH_HEADER": "secret",
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
        "AUTH_HEADER": "<redacted>",
        "AWS_ACCESS_KEY_ID": "<redacted>",
        "BATCH_KEY": "batch",
        "CREDENTIAL": "<redacted>",
        "GROUP_KEY": "sample",
        "LABEL_KEY": "cell_type",
        "OPENAI_API_KEY": "<redacted>",
        "OPENAI_KEY": "<redacted>",
        "PASSWORD_FILE": "<redacted>",
        "PRIVATE_KEY": "<redacted>",
        "VISIBLE": 1,
    }


def test_snapshot_effective_parameters_hides_absolute_paths(tmp_path: Path) -> None:
    project_root = tmp_path / "project"
    project_root.mkdir()

    assert snapshot_effective_parameters(
        {
            "INPUT_PATH": project_root / "data" / "input.h5ad",
            "RSCRIPT_BIN": str(tmp_path / "external" / "Rscript"),
        },
        path_root=project_root,
    ) == {
        "INPUT_PATH": "data/input.h5ad",
        "RSCRIPT_BIN": "<external>/Rscript",
    }


def test_snapshot_effective_parameters_rejects_unsupported_value() -> None:
    with pytest.raises(TypeError, match="UNSUPPORTED"):
        snapshot_effective_parameters({"UNSUPPORTED": {1, 2}})


def _init_git_repo(path: Path) -> Path:
    tracked = path / "tracked.txt"
    subprocess.run(["git", "init", "-q", str(path)], check=True)
    for key, value in (("user.email", "test@example.org"), ("user.name", "Test")):
        subprocess.run(["git", "-C", str(path), "config", key, value], check=True)
    tracked.write_text("before", encoding="utf-8")
    subprocess.run(["git", "-C", str(path), "add", "tracked.txt"], check=True)
    subprocess.run(["git", "-C", str(path), "commit", "-qm", "test"], check=True)
    return tracked


def test_collect_runtime_provenance_records_git_and_packages(tmp_path: Path) -> None:
    _init_git_repo(tmp_path)
    expected_commit = subprocess.check_output(
        ["git", "-C", str(tmp_path), "rev-parse", "HEAD"], text=True
    ).strip()

    provenance = collect_runtime_provenance(tmp_path, ("pytest", "definitely-missing-package"))

    assert provenance["python"] == sys.version.split()[0]
    assert provenance["platform"]
    assert provenance["git_commit"] == expected_commit
    assert provenance["git_available"] is True
    assert provenance["git_dirty"] is False
    assert provenance["git_status_sha256"] is None
    assert provenance["git_untracked_count"] == 0
    assert provenance["git_tracked_dirty"] is False
    assert provenance["git_tracked_diff_sha256"] is None
    assert provenance["packages"]["pytest"] != "unavailable"
    assert provenance["packages"]["definitely-missing-package"] == "unavailable"


def test_collect_runtime_provenance_marks_git_unavailable(tmp_path: Path) -> None:
    provenance = collect_runtime_provenance(tmp_path)
    assert provenance["git_available"] is False
    assert provenance["git_commit"] == "unavailable"
    unavailable = (
        "git_dirty", "git_status_sha256", "git_untracked_count",
        "git_tracked_dirty", "git_tracked_diff_sha256",
    )
    assert all(provenance[key] is None for key in unavailable)


def test_collect_runtime_provenance_hashes_tracked_diff(tmp_path: Path) -> None:
    tracked = _init_git_repo(tmp_path)
    tracked.write_text("after", encoding="utf-8")
    expected_diff = subprocess.check_output(
        ["git", "-C", str(tmp_path), "diff", "--binary", "HEAD"]
    )

    provenance = collect_runtime_provenance(tmp_path)
    assert provenance["git_dirty"] is True
    assert provenance["git_tracked_dirty"] is True
    assert provenance["git_tracked_diff_sha256"] == hashlib.sha256(expected_diff).hexdigest()


def test_collect_runtime_provenance_includes_untracked_files(tmp_path: Path) -> None:
    _init_git_repo(tmp_path)
    (tmp_path / "untracked.txt").write_text("new", encoding="utf-8")
    status = subprocess.check_output(
        ["git", "-C", str(tmp_path), "status", "--porcelain=v1", "--untracked-files=all", "-z"]
    )

    provenance = collect_runtime_provenance(tmp_path)
    assert provenance["git_dirty"] is True
    assert provenance["git_status_sha256"] == hashlib.sha256(status).hexdigest()
    assert provenance["git_untracked_count"] == 1
    assert provenance["git_tracked_dirty"] is False


def test_publish_compatibility_symlink_creates_and_reuses_relative_link(tmp_path: Path) -> None:
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

    assert publish_compatibility_symlink(link, target_one) == link
    with pytest.raises(FileExistsError):
        publish_compatibility_symlink(link, target_two)
    assert link.resolve() == target_one.resolve()


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


def test_publish_compatibility_symlink_does_not_overwrite_concurrent_entry(
    monkeypatch, tmp_path: Path
) -> None:
    target = tmp_path / "target.h5ad"
    target.write_text("target", encoding="utf-8")

    link = tmp_path / "latest.h5ad"

    def competing_create(self, target, *args, **kwargs):
        self.write_text("concurrent", encoding="utf-8")
        raise FileExistsError(self)

    monkeypatch.setattr(Path, "symlink_to", competing_create)
    with pytest.raises(FileExistsError):
        publish_compatibility_symlink(link, target)

    assert link.read_text(encoding="utf-8") == "concurrent"


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


def test_prepare_run_increment_preserves_digit_width(tmp_path: Path) -> None:
    """on_exists='increment' 时保持位数：run-001 → run-002 → run-003。"""
    p1 = prepare_run(tmp_path, "run-001", on_exists="increment")
    assert p1.run_id == "run-001"

    p2 = prepare_run(tmp_path, "run-001", on_exists="increment")
    assert p2.run_id == "run-002"

    p3 = prepare_run(tmp_path, "run-001", on_exists="increment")
    assert p3.run_id == "run-003"

    # 确认三个目录都存在
    assert (tmp_path / "run-001").is_dir()
    assert (tmp_path / "run-002").is_dir()
    assert (tmp_path / "run-003").is_dir()


def test_prepare_run_increment_appends_dash_for_no_digit_suffix(tmp_path: Path) -> None:
    """原 run_id 无数字后缀时追加 -2、-3。"""
    p1 = prepare_run(tmp_path, "my-run", on_exists="increment")
    assert p1.run_id == "my-run"

    p2 = prepare_run(tmp_path, "my-run", on_exists="increment")
    assert p2.run_id == "my-run-2"

    p3 = prepare_run(tmp_path, "my-run", on_exists="increment")
    assert p3.run_id == "my-run-3"

    assert (tmp_path / "my-run").is_dir()
    assert (tmp_path / "my-run-2").is_dir()
    assert (tmp_path / "my-run-3").is_dir()


def test_prepare_run_rejects_invalid_on_exists_value(tmp_path: Path) -> None:
    """非法 on_exists 值抛 ValueError，不静默降级。"""
    with pytest.raises(ValueError, match="on_exists must be"):
        prepare_run(tmp_path, "run-001", on_exists="overwrite")

    with pytest.raises(ValueError, match="on_exists must be"):
        prepare_run(tmp_path, "run-001", on_exists="")


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


def _artifact(paths, role: str, relative: str, content: bytes = b"artifact") -> dict:
    path = paths.draft_dir / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return {"role": role, "path": relative, "sha256": sha256_file(path)}


def _declare_artifacts(paths, artifacts) -> None:
    manifest = json.loads(paths.manifest_path.read_text(encoding="utf-8"))
    manifest["artifacts"] = artifacts
    atomic_write_json(paths.manifest_path, manifest, overwrite=True)


def test_legacy_manifest_without_artifacts_remains_valid(tmp_path: Path) -> None:
    paths, checkpoint = _write_draft(tmp_path, StageStatus.SUCCESS)
    assert validate_artifacts(paths.manifest_path) == {}
    assert validate_checkpoint(paths.manifest_path) == checkpoint
    assert resume_run(tmp_path, paths.run_id) == paths


def test_artifacts_validate_and_move_with_atomic_promotion(tmp_path: Path) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    _declare_artifacts(paths, [
        _artifact(paths, "report", "reports/report.json", b"report"),
        _artifact(paths, "table", "tables/table.csv", b"table"),
    ])
    validated = validate_artifacts(paths.manifest_path)
    assert validated == {
        "report": paths.draft_dir / "reports/report.json",
        "table": paths.draft_dir / "tables/table.csv",
    }
    assert resume_run(tmp_path, paths.run_id) == paths

    promote_run(paths)
    promoted_manifest = paths.promoted_dir / "manifest.json"
    assert validate_artifacts(promoted_manifest) == {
        "report": paths.promoted_dir / "reports/report.json",
        "table": paths.promoted_dir / "tables/table.csv",
    }
    assert resume_run(tmp_path, paths.run_id, promoted=True) == paths


@pytest.mark.parametrize("failure", ["tamper", "missing"])
def test_invalid_secondary_artifact_blocks_promotion_without_partial(
    tmp_path: Path, failure: str
) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    artifact = _artifact(paths, "report", "reports/report.json")
    _declare_artifacts(paths, [artifact])
    report = paths.draft_dir / artifact["path"]
    report.write_bytes(b"tampered") if failure == "tamper" else report.unlink()

    error = ValueError if failure == "tamper" else FileNotFoundError
    with pytest.raises(error):
        promote_run(paths)
    assert paths.draft_dir.is_dir() and not paths.promoted_dir.exists()


def test_promoted_resume_rejects_secondary_artifact_tampering(tmp_path: Path) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    artifact = _artifact(paths, "report", "reports/report.json")
    _declare_artifacts(paths, [artifact])
    promote_run(paths)
    (paths.promoted_dir / artifact["path"]).write_bytes(b"tampered")

    with pytest.raises(ValueError, match="artifact report hash"):
        resume_run(tmp_path, paths.run_id, promoted=True)


@pytest.mark.parametrize("value", [None, [], {}, "artifact"])
def test_artifacts_key_requires_nonempty_list(tmp_path: Path, value) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    _declare_artifacts(paths, value)
    with pytest.raises(ValueError, match="non-empty list"):
        validate_checkpoint(paths.manifest_path)


@pytest.mark.parametrize(("artifact", "match"), [
    ("not-an-object", "must be an object"),
    ({"path": "report.txt", "sha256": "0" * 64}, "role must"),
    ({"role": "report", "sha256": "0" * 64}, "path must"),
])
def test_artifact_item_requires_object_role_and_path(tmp_path: Path, artifact, match) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    _declare_artifacts(paths, [artifact])
    with pytest.raises(ValueError, match=match):
        validate_checkpoint(paths.manifest_path)


@pytest.mark.parametrize("duplicate", ["role", "path"])
def test_duplicate_artifact_roles_and_paths_are_rejected(
    tmp_path: Path, duplicate: str
) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    first = _artifact(paths, "report", "reports/report.json")
    second = _artifact(paths, "table", "tables/table.csv")
    if duplicate == "role":
        second["role"] = first["role"]
    else:
        second["path"], second["sha256"] = "reports//report.json", first["sha256"]
    _declare_artifacts(paths, [first, second])
    with pytest.raises(ValueError, match=f"duplicate artifact {duplicate}"):
        validate_checkpoint(paths.manifest_path)


@pytest.mark.parametrize("role", ["", "../report", "bad role", "/absolute"])
def test_artifact_role_must_be_nonempty_and_safe(tmp_path: Path, role: str) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    artifact = _artifact(paths, "report", "report.json")
    artifact["role"] = role
    _declare_artifacts(paths, [artifact])
    with pytest.raises(ValueError, match="role must be a non-empty safe name"):
        validate_checkpoint(paths.manifest_path)


@pytest.mark.parametrize("kind", ["traversal", "absolute"])
def test_artifact_path_cannot_escape_selected_run_state(tmp_path: Path, kind: str) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    outside = paths.run_dir / "outside.txt"
    outside.write_bytes(b"outside")
    raw_path = "../outside.txt" if kind == "traversal" else str(outside)
    _declare_artifacts(paths, [{
        "role": "report", "path": raw_path, "sha256": sha256_file(outside),
    }])
    with pytest.raises(ValueError, match="inside the run state directory"):
        validate_checkpoint(paths.manifest_path)


@pytest.mark.parametrize("raw_path", [".", "bad\nname.txt"])
def test_artifact_path_must_name_a_safe_file(tmp_path: Path, raw_path: str) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    _declare_artifacts(paths, [{"role": "report", "path": raw_path, "sha256": "0" * 64}])
    with pytest.raises(ValueError, match="safe relative file path"):
        validate_checkpoint(paths.manifest_path)


def test_artifact_symlink_is_rejected_even_when_target_hash_matches(tmp_path: Path) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    target = paths.draft_dir / "target.txt"
    target.write_bytes(b"target")
    link = paths.draft_dir / "report.txt"
    link.symlink_to(target.name)
    _declare_artifacts(paths, [{
        "role": "report", "path": link.name, "sha256": sha256_file(target),
    }])
    with pytest.raises(ValueError, match="contains a symlink"):
        validate_checkpoint(paths.manifest_path)


def test_artifact_directory_is_not_a_regular_file(tmp_path: Path) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    directory = paths.draft_dir / "report"
    directory.mkdir()
    _declare_artifacts(paths, [{"role": "report", "path": directory.name, "sha256": "0" * 64}])
    with pytest.raises(ValueError, match="must be a regular file"):
        validate_checkpoint(paths.manifest_path)


def test_manifest_symlink_cannot_redirect_resume_outside_run_state(tmp_path: Path) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    outside = tmp_path / "outside-manifest.json"
    outside.write_bytes(paths.manifest_path.read_bytes())
    paths.manifest_path.unlink()
    paths.manifest_path.symlink_to(outside)
    with pytest.raises(ValueError, match="manifest path must not be a symlink"):
        validate_checkpoint(paths.manifest_path)
    with pytest.raises(ValueError, match="manifest path must not be a symlink"):
        resume_run(tmp_path, paths.run_id)


def test_draft_symlink_cannot_redirect_validation_resume_or_promotion(
    tmp_path: Path,
) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    outside = tmp_path / "outside-draft"
    os.replace(paths.draft_dir, outside)
    sentinel = outside / "sentinel.txt"
    sentinel.write_bytes(b"unchanged")
    paths.draft_dir.symlink_to(outside, target_is_directory=True)

    for operation in (
        lambda: validate_checkpoint(paths.manifest_path),
        lambda: resume_run(tmp_path, paths.run_id),
        lambda: promote_run(paths),
    ):
        with pytest.raises(ValueError, match="run state path must not be a symlink"):
            operation()

    assert sentinel.read_bytes() == b"unchanged"
    assert paths.draft_dir.is_symlink()
    assert not paths.promoted_dir.exists()


def test_promoted_symlink_cannot_redirect_resume_or_promotion(tmp_path: Path) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    outside = tmp_path / "outside-promoted"
    outside.mkdir()
    sentinel = outside / "sentinel.txt"
    sentinel.write_bytes(b"unchanged")
    paths.promoted_dir.symlink_to(outside, target_is_directory=True)

    with pytest.raises(ValueError, match="run state path must not be a symlink"):
        resume_run(tmp_path, paths.run_id, promoted=True)
    with pytest.raises(ValueError, match="run state path must not be a symlink"):
        promote_run(paths)

    assert sentinel.read_bytes() == b"unchanged"
    assert paths.promoted_dir.is_symlink()
    assert paths.draft_dir.is_dir()


def test_run_root_symlink_is_rejected_before_resume(tmp_path: Path) -> None:
    real_root = tmp_path / "real-root"
    paths, _ = _write_draft(real_root, StageStatus.SUCCESS)
    linked_root = tmp_path / "linked-root"
    linked_root.symlink_to(real_root, target_is_directory=True)

    with pytest.raises(ValueError, match="run root path must not be a symlink"):
        resume_run(linked_root, paths.run_id)


def test_post_promotion_tampering_rolls_back_to_auditable_draft(
    monkeypatch, tmp_path: Path
) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    artifact = _artifact(paths, "report", "reports/report.json")
    _declare_artifacts(paths, [artifact])
    real_replace = os.replace

    def replace_then_tamper(source, destination):
        real_replace(source, destination)
        if Path(source) == paths.draft_dir.resolve():
            (Path(destination) / artifact["path"]).write_bytes(b"tampered")

    monkeypatch.setattr("scrna_integration.run_contract.os.replace", replace_then_tamper)

    with pytest.raises(ValueError, match="artifact report hash"):
        promote_run(paths)

    assert not paths.promoted_dir.exists()
    assert paths.draft_dir.is_dir()
    assert (paths.draft_dir / artifact["path"]).read_bytes() == b"tampered"
    with pytest.raises(ValueError, match="artifact report hash"):
        resume_run(tmp_path, paths.run_id)


def test_failed_post_promotion_rollback_reports_both_errors(
    monkeypatch, tmp_path: Path
) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    artifact = _artifact(paths, "report", "reports/report.json")
    _declare_artifacts(paths, [artifact])
    real_replace = os.replace
    calls = 0

    def replace_then_fail_rollback(source, destination):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("rollback denied")
        real_replace(source, destination)
        (Path(destination) / artifact["path"]).write_bytes(b"tampered")

    monkeypatch.setattr("scrna_integration.run_contract.os.replace", replace_then_fail_rollback)

    with pytest.raises(RuntimeError) as error:
        promote_run(paths)

    assert "artifact report hash does not match manifest" in str(error.value)
    assert "rollback denied" in str(error.value)


@pytest.mark.parametrize("sha256", [None, "bad", "0" * 64])
def test_artifact_hash_must_be_valid_and_live(tmp_path: Path, sha256) -> None:
    paths, _ = _write_draft(tmp_path, StageStatus.SUCCESS)
    artifact = _artifact(paths, "report", "report.txt")
    artifact["sha256"] = sha256
    _declare_artifacts(paths, [artifact])
    with pytest.raises(ValueError, match="sha256|hash does not match"):
        validate_checkpoint(paths.manifest_path)


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


# ---- validate_expression_contract ------------------------------------------------


class _FakeAdata:
    """最小 AnnData 替身，仅提供 uns dict 访问。"""

    def __init__(self, uns: dict | None = None) -> None:
        self.uns: dict = uns if uns is not None else {}


def _valid_contract(**overrides: object) -> dict:
    """返回一个通过所有 schema 校验的 expression_contract。"""
    contract: dict = {
        "x_scale": "raw_counts",
        "counts_layer": "counts",
        "counts_source": "X",
        "counts_validated": True,
        "counts_integer_check": "full",
        "soupx_layer": None,
        "processing_history": [],
        "stage": "01",
    }
    contract.update(overrides)
    return contract


def test_validate_expression_contract_valid_full() -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract()})
    result = validate_expression_contract(adata)
    assert result["x_scale"] == "raw_counts"
    assert result["counts_layer"] == "counts"


def test_validate_expression_contract_missing_key() -> None:
    contract = _valid_contract()
    del contract["counts_source"]
    adata = _FakeAdata({"expression_contract": contract})
    with pytest.raises(ValueError, match="missing required keys"):
        validate_expression_contract(adata)


def test_validate_expression_contract_unknown_key() -> None:
    contract = _valid_contract()
    contract["extra_field"] = "unexpected"
    adata = _FakeAdata({"expression_contract": contract})
    with pytest.raises(ValueError, match="has unknown keys"):
        validate_expression_contract(adata)


def test_validate_expression_contract_missing_entirely() -> None:
    adata = _FakeAdata({})
    with pytest.raises(KeyError, match="expression_contract not found"):
        validate_expression_contract(adata)


def test_validate_expression_contract_null_contract() -> None:
    adata = _FakeAdata({"expression_contract": None})
    with pytest.raises(KeyError, match="expression_contract not found"):
        validate_expression_contract(adata)


def test_validate_expression_contract_not_dict() -> None:
    adata = _FakeAdata({"expression_contract": "not a dict"})
    with pytest.raises(ValueError, match="must be a dict"):
        validate_expression_contract(adata)


@pytest.mark.parametrize("bad_scale", ["logcounts", "cpm", "raw", 1, None])
def test_validate_expression_contract_invalid_x_scale(bad_scale) -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(x_scale=bad_scale)})
    with pytest.raises(ValueError, match="x_scale"):
        validate_expression_contract(adata)


def test_validate_expression_contract_x_scale_matches_expected() -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(x_scale="raw_counts")})
    validate_expression_contract(adata, expected_scale="raw_counts")  # 不抛


def test_validate_expression_contract_x_scale_expected_mismatch() -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(x_scale="raw_counts")})
    with pytest.raises(ValueError, match="expected 'normalized_log1p'"):
        validate_expression_contract(adata, expected_scale="normalized_log1p")


@pytest.mark.parametrize("bad_source", [".raw", "layers.counts", "raw.X", "", 1, None])
def test_validate_expression_contract_invalid_counts_source(bad_source) -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(counts_source=bad_source)})
    with pytest.raises(ValueError, match="counts_source"):
        validate_expression_contract(adata)


@pytest.mark.parametrize("good_source", ["X", ".raw.X", "layers[counts]"])
def test_validate_expression_contract_valid_counts_sources(good_source) -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(counts_source=good_source)})
    validate_expression_contract(adata)  # 不抛


def test_validate_expression_contract_counts_validated_not_bool() -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(counts_validated="yes")})
    with pytest.raises(ValueError, match="counts_validated must be bool"):
        validate_expression_contract(adata)


@pytest.mark.parametrize("bad_check", ["partial", "sampled", None, 0])
def test_validate_expression_contract_invalid_integer_check(bad_check) -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(counts_integer_check=bad_check)})
    with pytest.raises(ValueError, match="counts_integer_check"):
        validate_expression_contract(adata)


def test_validate_expression_contract_valid_integer_checks() -> None:
    for check in ("full", "blockwise"):
        adata = _FakeAdata({"expression_contract": _valid_contract(counts_integer_check=check)})
        validate_expression_contract(adata)  # 不抛


@pytest.mark.parametrize("bad_soupx", ["counts_soupx_old", "soupx", 0, False])
def test_validate_expression_contract_invalid_soupx_layer(bad_soupx) -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(soupx_layer=bad_soupx)})
    with pytest.raises(ValueError, match="soupx_layer"):
        validate_expression_contract(adata)


def test_validate_expression_contract_soupx_layer_none_allowed() -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(soupx_layer=None)})
    validate_expression_contract(adata)  # 不抛


def test_validate_expression_contract_soupx_layer_counts_soupx_allowed() -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(soupx_layer="counts_soupx")})
    validate_expression_contract(adata)  # 不抛


def test_validate_expression_contract_processing_history_not_list() -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(processing_history="step1")})
    with pytest.raises(ValueError, match="processing_history must be a list"):
        validate_expression_contract(adata)


@pytest.mark.parametrize("bad_stage", ["04", "1", "A", "", None])
def test_validate_expression_contract_invalid_stage(bad_stage) -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(stage=bad_stage)})
    with pytest.raises(ValueError, match="stage must be one of"):
        validate_expression_contract(adata)


def test_validate_expression_contract_stage_matches_expected() -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(stage="01")})
    validate_expression_contract(adata, stage="01")  # 不抛


def test_validate_expression_contract_stage_expected_mismatch() -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(stage="02")})
    with pytest.raises(ValueError, match="expected '03'"):
        validate_expression_contract(adata, stage="03")


def test_validate_expression_contract_empty_counts_layer() -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(counts_layer="")})
    with pytest.raises(ValueError, match="counts_layer must be a non-empty string"):
        validate_expression_contract(adata)


def test_validate_expression_contract_returns_independent_dict() -> None:
    contract = _valid_contract()
    adata = _FakeAdata({"expression_contract": contract})
    result = validate_expression_contract(adata)
    assert result is not contract
    assert result == contract


def test_validate_expression_contract_x_scale_normalized_log1p() -> None:
    adata = _FakeAdata({"expression_contract": _valid_contract(x_scale="normalized_log1p")})
    validate_expression_contract(adata)

    with pytest.raises(ValueError, match="expected 'raw_counts'"):
        validate_expression_contract(adata, expected_scale="raw_counts")


# ---- aggregate_method_status ------------------------------------------------------


def test_aggregate_all_success() -> None:
    result = aggregate_method_status({
        "harmony": MethodStatus.SUCCESS,
        "scvi": MethodStatus.SUCCESS,
    })
    assert result == {"harmony": MethodStatus.SUCCESS, "scvi": MethodStatus.SUCCESS}


def test_aggregate_skipped_excluded() -> None:
    result = aggregate_method_status({
        "pca": MethodStatus.SUCCESS,
        "scvi": MethodStatus.SKIPPED_BY_USER,
        "harmony": MethodStatus.SUCCESS,
    })
    assert "scvi" not in result
    assert result["pca"] is MethodStatus.SUCCESS
    assert result["harmony"] is MethodStatus.SUCCESS


def test_aggregate_failed_included() -> None:
    result = aggregate_method_status({
        "scvi": MethodStatus.FAILED,
    })
    assert result == {"scvi": MethodStatus.FAILED}


def test_aggregate_unavailable_included() -> None:
    result = aggregate_method_status({
        "scanvi": MethodStatus.UNAVAILABLE,
    })
    assert result == {"scanvi": MethodStatus.UNAVAILABLE}


def test_aggregate_mixed_statuses() -> None:
    result = aggregate_method_status({
        "pca": MethodStatus.SUCCESS,
        "harmony": MethodStatus.SUCCESS,
        "scvi": MethodStatus.FAILED,
        "scanvi": MethodStatus.SKIPPED_BY_USER,
        "bbknn": MethodStatus.UNAVAILABLE,
    })
    assert result == {
        "pca": MethodStatus.SUCCESS,
        "harmony": MethodStatus.SUCCESS,
        "scvi": MethodStatus.FAILED,
        "bbknn": MethodStatus.UNAVAILABLE,
    }


def test_aggregate_empty_input() -> None:
    assert aggregate_method_status({}) == {}


def test_aggregate_all_skipped() -> None:
    result = aggregate_method_status({
        "scvi": MethodStatus.SKIPPED_BY_USER,
        "scanvi": MethodStatus.SKIPPED_BY_USER,
    })
    assert result == {}


def test_aggregate_accepts_string_status_values() -> None:
    result = aggregate_method_status({
        "pca": "success",
        "scvi": "failed",
        "harmony": "skipped_by_user",
    })
    assert result["pca"] is MethodStatus.SUCCESS
    assert result["scvi"] is MethodStatus.FAILED
    assert "harmony" not in result


def test_aggregate_result_usable_by_determine_stage_status() -> None:
    """aggregate_method_status 输出可直接传入 determine_stage_status。"""
    required = aggregate_method_status({
        "pca": MethodStatus.SUCCESS,
        "harmony": MethodStatus.SUCCESS,
        "scvi": MethodStatus.SKIPPED_BY_USER,
    })
    status = determine_stage_status(required, {"counts_valid": True})
    assert status is StageStatus.SUCCESS


def test_aggregate_with_failures_causes_determine_stage_status_failed() -> None:
    """含 FAILED 的聚合结果会让 determine_stage_status 返回 FAILED。"""
    required = aggregate_method_status({
        "pca": MethodStatus.SUCCESS,
        "scvi": MethodStatus.FAILED,
    })
    status = determine_stage_status(required, {"counts_valid": True})
    assert status is StageStatus.FAILED
