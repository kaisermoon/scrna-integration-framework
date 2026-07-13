"""PR1g1: versioned manifest v1 structural and state invariants."""

from __future__ import annotations

from pathlib import Path

import pytest

from scrna_integration.run_contract import (
    MANIFEST_SCHEMA_VERSION,
    MAX_MANIFEST_FAILURE_BYTES,
    atomic_write_json,
    validate_manifest,
)

MISSING = object()


def _manifest(status="SUCCESS"):
    manifest = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "run_id": "run-001",
        "stage": "04_embedded",
        "stage_status": status,
        "started_at": "2026-07-14T01:02:03Z",
        "completed_at": "2026-07-14T01:02:04.123Z",
        "inputs": [{"role": "upstream", "path": "upstream/manifest.json",
                    "kind": "file", "sha256": "1" * 64}],
        "effective_parameters": {},
        "runtime_provenance": {
            "python": "3.11", "platform": "test", "git_available": False, "git_commit": "unavailable",
            "git_dirty": None, "git_status_sha256": None, "git_untracked_count": None,
            "git_tracked_dirty": None, "git_tracked_diff_sha256": None, "packages": {},
        },
        "method_status": {"pca": "success"},
        "hard_postconditions": {"checkpoint_written": True},
        "warnings": ["warning"] if status == "SUCCESS_WITH_WARNINGS" else [],
        "artifacts": [],
        "checkpoint": {"path": "checkpoint.h5ad", "sha256": "0" * 64},
    }
    if status == "FAILED":
        manifest.pop("checkpoint")
        manifest["failure"] = {"type": "RuntimeError", "message": "failed"}
    return manifest


def _record(*, artifact=False, **updates):
    return {"role": "input", "path": "input/data.h5ad", "sha256": "1" * 64, **({} if artifact else {"kind": "file"})} | updates


def _available_runtime(**updates):
    return _manifest()["runtime_provenance"] | {"git_available": True, "git_commit": "a" * 40, "git_dirty": False, "git_status_sha256": None, "git_untracked_count": 0, "git_tracked_dirty": False, "git_tracked_diff_sha256": None} | updates


@pytest.mark.parametrize(("status", "state", "accepted"), [
    ("SUCCESS", "draft", False), ("SUCCESS", "promoted", False), ("SUCCESS_WITH_WARNINGS", "draft", False), ("SUCCESS_WITH_WARNINGS", "promoted", True), ("NEEDS_REVIEW", "draft", False), ("FAILED", "draft", False),
])
def test_valid_v1_states(status, state, accepted) -> None:
    manifest = _manifest(status)
    if accepted:
        manifest["warning_acceptance"] = {"accepted_by": "PI", "accepted_at": "2026-07-14T01:02:05Z"}
    if status == "SUCCESS" and state == "draft":
        manifest["runtime_provenance"] = _available_runtime(packages={"missing": "unavailable"})
    assert validate_manifest(manifest, state=state)["schema_version"] == "1"


@pytest.mark.parametrize(("key", "value"), [
    *[(key, MISSING) for key in ("schema_version", "run_id", "stage", "stage_status", "started_at", "completed_at", "inputs",
        "effective_parameters", "runtime_provenance", "method_status", "hard_postconditions", "warnings", "artifacts", "checkpoint")],
    ("schema_version", 1), ("schema_version", "2"), ("run_id", "bad id"), ("run_id", "运行"),
    ("stage", "bad/stage"), ("inputs", []), ("inputs", [{}]), ("inputs", ["input"]), ("effective_parameters", []),
    ("runtime_provenance", []), ("method_status", []), ("hard_postconditions", []), ("warnings", {}), ("warnings", [1]),
    ("artifacts", {}), ("artifacts", [{}]), ("artifacts", ["artifact"]), ("checkpoint", []), ("checkpoint", {}),
    ("effective_parameters", {"x": float("nan")}), ("effective_parameters", {"x": float("inf")}), ("effective_parameters", {"x": Path("x")}),
    ("effective_parameters", {"x": {1}}), ("effective_parameters", {1: "x"}), ("runtime_provenance", {}),
    ("runtime_provenance", _available_runtime(git_commit="bad")), ("runtime_provenance", _available_runtime(git_available=1)), ("runtime_provenance", _available_runtime(git_dirty=1)), ("runtime_provenance", _available_runtime(git_untracked_count=True)),
    ("runtime_provenance", _available_runtime(git_untracked_count=-1)), ("runtime_provenance", _available_runtime(git_status_sha256="0" * 64)), ("runtime_provenance", _available_runtime(git_dirty=True, git_status_sha256="bad")),
    ("runtime_provenance", _available_runtime(packages={1: "version"})), ("runtime_provenance", _available_runtime(packages={"pkg": 1})), ("runtime_provenance", _available_runtime(packages={"pkg": None})), ("runtime_provenance", dict(_available_runtime(), unknown="x")),
    ("method_status", {}), ("method_status", {"pca": 1}), ("method_status", {1: "success"}), ("method_status", {"": "success"}), ("method_status", {"pca": ""}), ("method_status", {"pca": "banana"}),
    ("hard_postconditions", {}), ("hard_postconditions", {"gate": 1}), ("hard_postconditions", {1: True}), ("warnings", [""]), ("warnings", ["same", "same"]),
    ("inputs", [_record(role="bad role")]), ("inputs", [_record(path="../escape")]), ("inputs", [_record(kind="blob")]), ("inputs", [_record(sha256="bad")]),
    ("inputs", [_record(path="a"), _record(path="b")]), ("inputs", [_record(), _record(role="other")]),
    ("artifacts", [_record(artifact=True, role=None)]), ("artifacts", [_record(artifact=True, path="/absolute")]), ("artifacts", [_record(artifact=True, sha256="bad")]),
    ("artifacts", [_record(artifact=True, path="checkpoint.h5ad")]), ("artifacts", [_record(artifact=True, path="a"), _record(artifact=True, path="b")]),
    ("artifacts", [_record(artifact=True), _record(artifact=True, role="other")]),
    ("checkpoint", {"path": "../escape", "sha256": "0" * 64}), ("checkpoint", {"path": "checkpoint.h5ad", "sha256": "bad"}),
    ("started_at", None), ("started_at", "2026-07-14 01:02:03Z"), ("started_at", "2026-02-30T01:02:03Z"),
    ("started_at", "2026-07-14T01:02:03+00:00"), ("completed_at", "2026-07-14T01:02:02Z"), ("stage_status", "UNKNOWN"), ("unknown_top_level", 1),
    ("checkpont", {"path": "x", "sha256": "0" * 64}), ("failure", {"message": "forbidden"}),
])
def test_field_and_record_container_types(key, value) -> None:
    manifest = _manifest()
    if value is MISSING:
        manifest.pop(key)
    else:
        manifest[key] = value
    with pytest.raises(ValueError):
        validate_manifest(manifest, state="draft")


@pytest.mark.parametrize(("status", "state", "warning_values", "acceptance"), [
    ("SUCCESS", "draft", ["unexpected"], MISSING), ("SUCCESS", "draft", [], {"accepted_by": "PI", "accepted_at": "2026-07-14T01:02:05Z"}),
    ("SUCCESS_WITH_WARNINGS", "draft", [], MISSING), ("SUCCESS_WITH_WARNINGS", "draft", ["warning"], {"accepted_by": ""}),
    ("NEEDS_REVIEW", "draft", [], {"accepted_by": "PI", "accepted_at": "2026-07-14T01:02:05Z"}),
    ("FAILED", "draft", [], {"accepted_by": "PI", "accepted_at": "2026-07-14T01:02:05Z"}),
    *[("SUCCESS_WITH_WARNINGS", "promoted", ["warning"], acceptance) for acceptance in
      (MISSING, {}, {"accepted_by": "PI"}, "accepted", {"accepted_by": "PI", "accepted_at": "2026-07-14"},
       {"accepted_by": "PI", "accepted_at": "2026-07-14T01:02:02Z"})],
    ("FAILED", "promoted", [], MISSING), ("NEEDS_REVIEW", "promoted", [], MISSING),
])
def test_status_warning_semantics(status, state, warning_values, acceptance) -> None:
    manifest = _manifest(status)
    manifest["warnings"] = warning_values
    if acceptance is not MISSING:
        manifest["warning_acceptance"] = acceptance
    with pytest.raises(ValueError):
        validate_manifest(manifest, state=state)


@pytest.mark.parametrize("failure", [MISSING, {}, [], {"message": float("nan")}])
def test_failed_manifest_requires_finite_bounded_failure(failure) -> None:
    manifest = _manifest("FAILED")
    if failure is MISSING:
        manifest.pop("failure")
    else:
        manifest["failure"] = failure
    with pytest.raises(ValueError):
        validate_manifest(manifest, state="draft")


def test_failed_manifest_forbids_checkpoint_but_allows_diagnostics() -> None:
    manifest = _manifest("FAILED")
    manifest["checkpoint"] = {"path": "partial.h5ad"}
    with pytest.raises(ValueError, match="status-forbidden"):
        validate_manifest(manifest, state="draft")
    manifest.pop("checkpoint")
    manifest["failure"] = {"message": "x" * MAX_MANIFEST_FAILURE_BYTES}
    with pytest.raises(ValueError, match="size limit"):
        validate_manifest(manifest, state="draft")
    manifest["failure"] = {"message": "failed"}
    manifest["artifacts"] = [_record(artifact=True, role="diagnostic", path="diagnostic/report.json")]
    assert validate_manifest(manifest, state="draft")["artifacts"]


def test_legacy_requires_explicit_compatibility_path() -> None:
    legacy = {"run_id": "legacy", "stage_status": "SUCCESS"}
    with pytest.raises(ValueError, match="state is required"):
        validate_manifest(_manifest())
    with pytest.raises(ValueError, match="draft or promoted"):
        validate_manifest(_manifest(), state="review")
    with pytest.raises(ValueError, match="schema_version"):
        validate_manifest(legacy, state="draft")
    with pytest.warns(DeprecationWarning):
        assert validate_manifest(legacy, state="draft", allow_legacy=True) == legacy


def test_path_context_is_strict(tmp_path: Path) -> None:
    path = tmp_path / "run-001" / "draft" / "manifest.json"
    path.parent.mkdir(parents=True)
    atomic_write_json(path, _manifest())
    assert validate_manifest(path)["run_id"] == "run-001"
    with pytest.raises(ValueError, match="state does not match"):
        validate_manifest(path, state="promoted")
    with pytest.raises(ValueError, match="expected run_id"):
        validate_manifest(path, expected_run_id="other")


def test_runtime_provenance_accepts_strict_r_environment() -> None:
    manifest = _manifest()
    manifest["runtime_provenance"]["r_environment"] = {
        "available": True, "version": "4.4.1", "packages": {"Seurat": "5.1.0"},
    }
    assert validate_manifest(manifest, state="draft")["runtime_provenance"]["r_environment"]["available"]
    manifest["runtime_provenance"]["r_environment"]["error"] = "forbidden"
    with pytest.raises(ValueError, match="r_environment"):
        validate_manifest(manifest, state="draft")
