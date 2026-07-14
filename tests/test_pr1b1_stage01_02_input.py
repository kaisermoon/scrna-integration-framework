from __future__ import annotations

import ast
import gc
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import scrna_integration.run_contract as rc
ROOT = Path(__file__).parents[1]
SOURCES = {
    "01_kim.ipynb": "Kim_2023", "01_nancang.ipynb": "Nancang_2025",
    "01_nowicki.ipynb": "Nowicki_2023", "01_yue.ipynb": "Yue_2024",
}
NB01 = [f"notebooks/01_per_dataset/{name}" for name in SOURCES]
NB02 = "notebooks/02_merged.ipynb"
def _nb(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))

def _source(cell: dict) -> str:
    source = cell.get("source", "")
    return "".join(source) if isinstance(source, list) else source

def _cell(path: str, marker: str) -> str:
    for cell in _nb(path)["cells"]:
        if marker in (source := _source(cell)):
            return source
    raise AssertionError(f"{path} 缺少 cell: {marker}")

class _Adata:
    def __init__(self, source: str, n_obs: int = 2) -> None:
        self.X = sp.csr_matrix(np.ones((n_obs, 2), dtype=np.float32))
        self.n_obs, self.n_vars = self.X.shape
        self.var_names = ["GENE1", "GENE2"]
        self.obs = pd.DataFrame({"source_dataset": [source] * n_obs})
        self.uns: dict = {}

    def write_h5ad(self, path: Path, **_: object) -> None:
        Path(path).write_bytes(b"checkpoint")

    def var_names_make_unique(self) -> None:
        pass

def _env01(tmp: Path, path: str, run_id: str, n_obs: int = 2) -> dict:
    env: dict = {}
    exec(_cell(path, "# === PARAMS ==="), env)
    source = SOURCES[Path(path).name]
    manifest = tmp / f"{source}.yaml"
    manifest.write_text(f"source_dataset: {source}\npreprocessing_done: []\n", encoding="utf-8")
    env.update(
        adata=_Adata(source, n_obs), source_dataset=source, MANIFEST_PATH=str(manifest),
        RUN_ID=run_id, RUN_ROOT=str(tmp / "runs"), OUTPUT_FILENAME=f"{source}.h5ad",
        _root=str(tmp), os=os, sp=sp, np=np, Path=Path, gc=gc,
    )
    for name in (
        "atomic_write_json", "collect_runtime_provenance", "determine_stage_status",
        "prepare_run", "promote_run", "sha256_file", "snapshot_effective_parameters",
    ):
        env[name] = getattr(rc, name)
    return env

def _upstream(root: Path, source: str, run_id: str) -> Path:
    paths = rc.prepare_run(root, run_id)
    checkpoint = paths.draft_dir / f"{source}.h5ad"
    checkpoint.write_bytes(source.encode())
    rc.atomic_write_json(paths.manifest_path, {
        "run_id": run_id, "stage": "01_qcd", "stage_status": "SUCCESS",
        "source_dataset": source,
        "checkpoint": {"path": checkpoint.name, "sha256": rc.sha256_file(checkpoint)},
    })
    return rc.promote_run(paths)

def _env02(tmp: Path) -> tuple[dict, dict[Path, _Adata]]:
    root = tmp / "runs"
    runs = {source: f"run-{i}" for i, source in enumerate(SOURCES.values())}
    registry = {_upstream(root, src, rid).resolve(): _Adata(src) for src, rid in runs.items()}
    return {
        "UPSTREAM_RUN_ROOT": str(root), "UPSTREAM_RUNS": runs, "Path": Path, "json": json,
        "pd": pd, "display": lambda _: None, "resume_run": rc.resume_run,
        "sha256_file": rc.sha256_file, "validate_checkpoint": rc.validate_checkpoint,
        "sc": SimpleNamespace(read_h5ad=lambda path: registry[Path(path).resolve()]),
    }, registry

def test_json_ast_params_and_prepare_run_placement() -> None:
    for path in [*NB01, NB02]:
        code = [_source(c) for c in _nb(path)["cells"] if c["cell_type"] == "code"]
        ast.parse("\n".join(code), filename=path)
    for path in NB01:
        assert all(name in _cell(path, "# === PARAMS ===") for name in ("RUN_ID", "RUN_ROOT", "OUTPUT_FILENAME"))
        final = _cell(path, "# Checkpoint")
        contract_cell = _cell(path, "# === 数据读入") if path == NB01[0] else final
        assert "snapshot_effective_parameters(globals()" in contract_cell and "collect_runtime_provenance" in contract_cell
        earlier = [_source(c) for c in _nb(path)["cells"] if c["cell_type"] == "code" and _source(c) != final]
        if path != NB01[0]:
            assert not any(isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "prepare_run"
                           for source in earlier for n in ast.walk(ast.parse(source)))
    assert all(" in dir()" not in _cell(path, "qc_report = {") for path in (NB01[0], NB01[3]))

@pytest.mark.parametrize("path", NB01[1:])
def test_stage01_success_promotes_auditable_checkpoint(tmp_path: Path, path: str) -> None:
    env = _env01(tmp_path, path, f"success-{Path(path).stem}")
    expected = env["source_dataset"]
    exec(_cell(path, "# Checkpoint"), env)
    promoted = tmp_path / "runs" / env["RUN_ID"] / "promoted"
    manifest = json.loads((promoted / "manifest.json").read_text(encoding="utf-8"))
    assert (manifest["stage"], manifest["source_dataset"], manifest["stage_status"]) == ("01_qcd", expected, "SUCCESS")
    assert manifest["effective_parameters"]["RUN_ID"] == env["RUN_ID"]
    assert not {"RSCRIPT_BIN", "R_AVAILABLE"} & manifest["effective_parameters"].keys()
    assert manifest["runtime_provenance"]["python"] and all(manifest["hard_postconditions"].values())
    assert rc.validate_checkpoint(promoted / "manifest.json").is_file()

@pytest.mark.parametrize("values", [[], ["Wrong", "Wrong"], ["Nancang_2025", "Other"], ["Nancang_2025", np.nan]])
def test_stage01_bad_source_or_empty_fails_without_checkpoint(tmp_path: Path, values: list) -> None:
    env = _env01(tmp_path, NB01[1], "failed", n_obs=len(values))
    env["adata"].obs["source_dataset"] = values
    with pytest.raises(RuntimeError, match="Stage 01 FAILED"):
            exec(_cell(NB01[1], "# Checkpoint"), env)
    draft = tmp_path / "runs/failed/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage_status"] == "FAILED" and "checkpoint" not in manifest
    assert not list(draft.glob("*.h5ad"))

def test_stage01_missing_manifest_does_not_claim_run_id(tmp_path: Path) -> None:
    env = _env01(tmp_path, NB01[1], "missing")
    Path(env["MANIFEST_PATH"]).unlink()
    with pytest.raises(FileNotFoundError):
        exec(_cell(NB01[1], "# Checkpoint"), env)
    assert not Path(env["RUN_ROOT"]).exists()

def test_stage01_partial_write_is_removed_and_audited(tmp_path: Path) -> None:
    env = _env01(tmp_path, NB01[1], "write-failed")
    def fail_write(path: Path, **_: object) -> None:
        Path(path).write_bytes(b"partial")
        raise OSError("disk full")
    env["adata"].write_h5ad = fail_write
    with pytest.raises(OSError, match="disk full"):
        exec(_cell(NB01[1], "# Checkpoint"), env)
    draft = tmp_path / "runs/write-failed/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage_status"] == "FAILED" and manifest["failure"] == {"type": "OSError", "message": "disk full"}
    assert "checkpoint" not in manifest and not list(draft.glob("*.h5ad"))

def test_stage02_loads_four_verified_sources(tmp_path: Path) -> None:
    env, _ = _env02(tmp_path)
    exec(_cell(NB02, "# 加载每个 per-dataset h5ad"), env)
    assert len(env["adatas"]) == 4
    assert {item["source_dataset"] for item in env["upstream_inputs"]} == set(SOURCES.values())
    assert all(item["manifest_sha256"] and item["checkpoint_sha256"] for item in env["upstream_inputs"])

@pytest.mark.parametrize(("case", "match"), [
    ("duplicate", "run ID.*唯一"), ("stage", "stage"), ("hash", "hash"), ("obs", "obs source_dataset"),
])
def test_stage02_rejects_invalid_upstream(tmp_path: Path, case: str, match: str) -> None:
    env, registry = _env02(tmp_path)
    if case == "duplicate":
        env["UPSTREAM_RUNS"]["Kim_2023"] = env["UPSTREAM_RUNS"]["Nancang_2025"]
    elif case == "stage":
        path = Path(env["UPSTREAM_RUN_ROOT"]) / env["UPSTREAM_RUNS"]["Kim_2023"] / "promoted/manifest.json"
        manifest = json.loads(path.read_text(encoding="utf-8")); manifest["stage"] = "02_merged"
        rc.atomic_write_json(path, manifest, overwrite=True)
    elif case == "hash":
        next(iter(registry)).write_bytes(b"tampered")
    else:
        path = next(path for path in registry if path.read_bytes() == b"Kim_2023")
        registry[path] = _Adata("Wrong_source")
    with pytest.raises(ValueError, match=match):
        exec(_cell(NB02, "# 加载每个 per-dataset h5ad"), env)
