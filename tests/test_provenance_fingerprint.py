from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pytest

import scrna_integration.run_contract as contract


def test_utc_now_is_rfc3339_utc() -> None:
    value = contract.utc_now_rfc3339()
    assert value.endswith("Z") and contract._RFC3339_UTC_RE.fullmatch(value)


def test_file_fingerprint_is_deterministic(tmp_path: Path) -> None:
    source = tmp_path / "input.bin"
    source.write_bytes(b"counts")
    expected = hashlib.sha256(b"counts").hexdigest()
    record = contract.fingerprint_input("counts", source, tmp_path)
    assert record == {"role": "counts", "path": "input.bin", "kind": "file", "sha256": expected}
    assert contract.fingerprint_input("counts", source, tmp_path)["sha256"] == expected


def test_tree_hash_uses_sorted_paths(tmp_path: Path) -> None:
    tree = tmp_path / "tree"
    tree.mkdir()
    (tree / "z.txt").write_text("z", encoding="utf-8")
    (tree / "a.txt").write_text("a", encoding="utf-8")
    records = [(name, 1, hashlib.sha256(value).hexdigest())
               for name, value in (("a.txt", b"a"), ("z.txt", b"z"))]
    canonical = json.dumps(records, ensure_ascii=False, separators=(",", ":")).encode()
    expected = hashlib.sha256(canonical).hexdigest()
    assert contract.fingerprint_input("matrix", tree, tmp_path, kind="tree")["sha256"] == expected


@pytest.mark.parametrize("mutation", ["bytes", "path", "size", "empty_directory"])
def test_tree_detects_mutation_between_scans(monkeypatch, tmp_path: Path, mutation: str) -> None:
    tree = tmp_path / "tree"
    tree.mkdir()
    source = tree / "a"
    source.write_bytes(b"one")
    original = contract._scan_tree
    calls = 0

    def mutate(path):
        nonlocal calls
        result = original(path)
        calls += 1
        if calls == 1:
            if mutation == "bytes":
                source.write_bytes(b"two")
            elif mutation == "path":
                source.rename(tree / "b")
            elif mutation == "size":
                source.write_bytes(b"longer")
            else:
                (tree / "empty").mkdir()
        return result

    monkeypatch.setattr(contract, "_scan_tree", mutate)
    with pytest.raises(RuntimeError, match="tree changed"):
        contract.fingerprint_input("tree", tree, tmp_path, kind="tree")


def test_rejects_symlinks_and_root_escape(tmp_path: Path) -> None:
    root = tmp_path / "root"
    root.mkdir()
    outside = tmp_path / "outside"
    outside.write_text("outside", encoding="utf-8")
    link = root / "link"
    link.symlink_to(outside)
    with pytest.raises(ValueError, match="regular non-symlink"):
        contract.fingerprint_input("link", link, root)
    with pytest.raises(ValueError, match="inside root"):
        contract.fingerprint_input("escape", outside, root)


def test_tree_rejects_case_collisions(tmp_path: Path) -> None:
    tree = tmp_path / "tree"
    tree.mkdir()
    (tree / "A").write_text("upper", encoding="utf-8")
    lower = tree / "a"
    lower.write_text("lower", encoding="utf-8")
    if len(list(tree.iterdir())) < 2:
        pytest.skip("filesystem is case-insensitive")
    with pytest.raises(ValueError, match="case-colliding"):
        contract.fingerprint_input("tree", tree, tmp_path, kind="tree")


def test_file_fingerprint_detects_toctou(monkeypatch, tmp_path: Path) -> None:
    source = tmp_path / "input"
    source.write_text("stable", encoding="utf-8")
    original = contract._stat_identity
    calls = 0

    def changing(info):
        nonlocal calls
        calls += 1
        return (*original(info), calls)

    monkeypatch.setattr(contract, "_stat_identity", changing)
    with pytest.raises(RuntimeError, match="changed"):
        contract.fingerprint_input("input", source, tmp_path)


def test_artifact_record_is_confined(tmp_path: Path) -> None:
    state = tmp_path / "draft"
    state.mkdir()
    artifact = state / "plot.pdf"
    artifact.write_bytes(b"plot")
    assert contract.artifact_record("plot", artifact, state)["path"] == "plot.pdf"
    with pytest.raises(ValueError, match="inside root"):
        contract.artifact_record("outside", tmp_path / "outside", state)


def test_collect_r_environment_success(monkeypatch) -> None:
    def run(command, **kwargs):
        assert isinstance(command, list) and kwargs["timeout"] == 3
        return subprocess.CompletedProcess(command, 0, "R\t4.4.1\nP\tSeurat\t5.1.0\n", "")

    monkeypatch.setattr(contract.subprocess, "run", run)
    expected = {"available": True, "version": "4.4.1", "packages": {"Seurat": "5.1.0"}}
    assert contract.collect_r_environment(["Seurat"], timeout=3) == expected


@pytest.mark.parametrize(
    ("error", "message"),
    [(FileNotFoundError(), "unavailable"), (PermissionError(), "unavailable"), (subprocess.TimeoutExpired("hidden", 1), "timed out")],
)
def test_collect_r_environment_unavailable(monkeypatch, error, message) -> None:
    monkeypatch.setattr(contract.subprocess, "run", lambda *args, **kwargs: (_ for _ in ()).throw(error))
    result = contract.collect_r_environment(["Seurat"], rscript="/secret/Rscript")
    assert result["available"] is False and message in result["error"] and "/secret" not in result["error"]


def test_collect_r_environment_rejects_malformed_output(monkeypatch) -> None:
    monkeypatch.setattr(contract.subprocess, "run", lambda *args, **kwargs:
                        subprocess.CompletedProcess([], 0, "not provenance\n", "credential=secret"))
    environment = contract.collect_r_environment(["Seurat"])
    assert environment["available"] is False and environment["error"] == "malformed R provenance output"
