from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pytest

import scrna_integration.run_contract as contract


class TestProvenanceFingerprint:
    def test_utc_now_is_rfc3339_utc(self) -> None:
        value = contract.utc_now_rfc3339()
        assert value.endswith("Z") and contract._RFC3339_UTC_RE.fullmatch(value)

    def test_file_fingerprint_is_deterministic(self, tmp_path: Path) -> None:
        (source := tmp_path / "input.bin").write_bytes(b"counts")
        expected = hashlib.sha256(b"counts").hexdigest()
        record = contract.fingerprint_input("counts", source, tmp_path)
        assert record == {"role": "counts", "path": "input.bin", "kind": "file", "sha256": expected}
        assert contract.fingerprint_input("counts", source, tmp_path)["sha256"] == expected

    def test_tree_hash_uses_sorted_paths(self, tmp_path: Path) -> None:
        (tree := tmp_path / "tree").mkdir()
        for name, value in (("z.txt", b"z"), ("a.txt", b"a")):
            (tree / name).write_bytes(value)
        records = [("file", name, 1, hashlib.sha256(value).hexdigest())
                   for name, value in (("a.txt", b"a"), ("z.txt", b"z"))]
        canonical = json.dumps(records, ensure_ascii=False, separators=(",", ":")).encode()
        assert contract.fingerprint_input("matrix", tree, tmp_path, kind="tree")["sha256"] == hashlib.sha256(canonical).hexdigest()

    @pytest.mark.parametrize("action", ["add", "remove", "rename"])
    def test_empty_directories_change_stable_tree_digest(self, tmp_path: Path, action: str) -> None:
        (tree := tmp_path / "tree").mkdir()
        empty = tree / "empty"
        empty.mkdir() if action != "add" else None
        before = contract.fingerprint_input("tree", tree, tmp_path, kind="tree")["sha256"]
        empty.mkdir() if action == "add" else empty.rmdir() if action == "remove" else empty.rename(tree / "renamed")
        assert contract.fingerprint_input("tree", tree, tmp_path, kind="tree")["sha256"] != before

    @pytest.mark.parametrize("mutation", ["bytes", "path", "size", "empty_directory"])
    def test_tree_detects_mutation_between_scans(self, monkeypatch, tmp_path: Path, mutation: str) -> None:
        (tree := tmp_path / "tree").mkdir()
        (source := tree / "a").write_bytes(b"one")
        actions = {"bytes": lambda: source.write_bytes(b"two"),
                   "path": lambda: source.rename(tree / "b"),
                   "size": lambda: source.write_bytes(b"longer"),
                   "empty_directory": lambda: (tree / "empty").mkdir()}
        original, calls = contract._scan_tree, 0

        def mutate(path):
            nonlocal calls
            result = original(path)
            if calls == 0:
                actions[mutation]()
            calls += 1
            return result
        monkeypatch.setattr(contract, "_scan_tree", mutate)
        with pytest.raises(RuntimeError, match="tree changed"):
            contract.fingerprint_input("tree", tree, tmp_path, kind="tree")

    def test_rejects_symlinks_and_root_escape(self, tmp_path: Path) -> None:
        (root := tmp_path / "root").mkdir()
        (outside := tmp_path / "outside").write_text("outside", encoding="utf-8")
        (link := root / "link").symlink_to(outside)
        with pytest.raises(ValueError, match="regular non-symlink"):
            contract.fingerprint_input("link", link, root)
        with pytest.raises(ValueError, match="inside root"):
            contract.fingerprint_input("escape", outside, root)

    def test_tree_rejects_case_collisions(self, tmp_path: Path) -> None:
        (tree := tmp_path / "tree").mkdir()
        (tree / "A").write_text("upper", encoding="utf-8")
        (tree / "a").write_text("lower", encoding="utf-8")
        if len(list(tree.iterdir())) < 2:
            pytest.skip("filesystem is case-insensitive")
        with pytest.raises(ValueError, match="case-colliding"):
            contract.fingerprint_input("tree", tree, tmp_path, kind="tree")

    def test_file_fingerprint_detects_toctou(self, monkeypatch, tmp_path: Path) -> None:
        (source := tmp_path / "input").write_text("stable", encoding="utf-8")
        original, identities = contract._stat_identity, iter(range(20))
        monkeypatch.setattr(contract, "_stat_identity", lambda info: (*original(info), next(identities)))
        with pytest.raises(RuntimeError, match="changed"):
            contract.fingerprint_input("input", source, tmp_path)

    def test_artifact_record_is_confined(self, tmp_path: Path) -> None:
        (state := tmp_path / "draft").mkdir()
        (artifact := state / "plot.pdf").write_bytes(b"plot")
        assert contract.artifact_record("plot", artifact, state)["size"] == 4
        with pytest.raises(ValueError, match="inside root"):
            contract.artifact_record("outside", tmp_path / "outside", state)

    def test_collect_r_environment_success(self, monkeypatch) -> None:
        def run(command, **kwargs):
            assert isinstance(command, list) and kwargs["timeout"] == 3
            return subprocess.CompletedProcess(command, 0, "R\t4.4.1\nP\tSeurat\t5.1.0\n", "")
        monkeypatch.setattr(contract.subprocess, "run", run)
        expected = {"available": True, "version": "4.4.1", "packages": {"Seurat": "5.1.0"}}
        assert contract.collect_r_environment(["Seurat"], timeout=3) == expected

    @pytest.mark.parametrize(("error", "message"), [(FileNotFoundError(), "unavailable"),
                             (PermissionError(), "unavailable"),
                             (subprocess.TimeoutExpired("hidden", 1), "timed out")])
    def test_collect_r_environment_unavailable(self, monkeypatch, error, message) -> None:
        monkeypatch.setattr(contract.subprocess, "run", lambda *args, **kwargs: (_ for _ in ()).throw(error))
        result = contract.collect_r_environment(["Seurat"], rscript="/secret/Rscript")
        assert result["available"] is False and message in result["error"] and "/secret" not in result["error"]

    def test_collect_r_environment_rejects_malformed_output(self, monkeypatch) -> None:
        monkeypatch.setattr(contract.subprocess, "run", lambda *args, **kwargs:
                            subprocess.CompletedProcess([], 0, "not provenance\n", "credential=secret"))
        environment = contract.collect_r_environment(["Seurat"])
        assert environment["available"] is False and environment["error"] == "malformed R provenance output"
