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


def _params_source(path: str) -> str:
    """返回 notebook 中全部 PARAMS 参数定义的拼接源码（兼容新旧 PARAMS 结构）。

    P1-e 将 01_nancang 的单 cell PARAMS 拆为四组（每组 1 md header + 1 code cell）；
    其他 notebook 仍使用旧结构（单 # === PARAMS === cell）。此函数自动适配两种结构。
    """
    nb = _nb(path)
    group_markers = [
        "### 1. 数据源", "### 2. QC 阈值",
        "### 3. 方法开关", "### 4. 输出版本与运行标识",
    ]
    found_group = False
    sources = []
    for i, cell in enumerate(nb["cells"]):
        src = _source(cell)
        if any(m in src for m in group_markers) and cell.get("cell_type") == "markdown":
            found_group = True
            # 下一个 cell 是对应的 code cell
            if i + 1 < len(nb["cells"]) and nb["cells"][i + 1].get("cell_type") == "code":
                sources.append(_source(nb["cells"][i + 1]))
    if found_group and sources:
        return "\n".join(sources)
    # fallback: 旧结构单 # === PARAMS === cell
    return _cell(path, "# === PARAMS ===")

class _Adata:
    def __init__(self, source: str, n_obs: int = 2, layers: dict | None = None) -> None:
        self.X = sp.csr_matrix(np.ones((n_obs, 2), dtype=np.float32))
        self.n_obs, self.n_vars = self.X.shape
        self.shape = self.X.shape
        self.var_names = ["GENE1", "GENE2"]
        self.obs = pd.DataFrame({"source_dataset": [source] * n_obs})
        self.uns: dict = {}
        self.layers: dict = {}

    def write_h5ad(self, path: Path, **_: object) -> None:
        Path(path).write_bytes(b"checkpoint")

    def var_names_make_unique(self) -> None:
        pass

def _env01(tmp: Path, path: str, run_id: str, n_obs: int = 2) -> dict:
    env: dict = {}
    exec(_params_source(path), env)
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
        "validate_expression_contract",
    ):
        env[name] = getattr(rc, name)
    # 填充 layers["counts"] 与 expression_contract，供 checkpoint 断言消费
    # 按 source 区分契约值：nowicki X 是 normalized_log1p 且 counts 来自 .raw.X
    # 四个 notebook 的 checkpoint 均调用 validate_expression_contract(full schema)，
    # 因此 counts_integer_check 不能为 None，必须为 "full" 或 "blockwise"
    _a = env["adata"]
    _a.layers["counts"] = sp.csr_matrix(np.ones((n_obs, 2), dtype=np.float32))
    _is_nowicki = source == "Nowicki_2023"
    _a.uns["expression_contract"] = {
        "x_scale": "normalized_log1p" if _is_nowicki else "raw_counts",
        "counts_layer": "counts",
        "counts_source": ".raw.X" if _is_nowicki else "X",
        "counts_validated": False,
        "counts_integer_check": "blockwise" if _is_nowicki else "full",
        "soupx_layer": None,
        "processing_history": [],
        "stage": "01",
    }
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
    registry = {}
    for src, rid in runs.items():
        ad = _Adata(src)
        ad.layers["counts"] = sp.csr_matrix(np.ones((ad.n_obs, 2), dtype=np.float32))
        ad.uns["expression_contract"] = {
            "x_scale": "raw_counts",
            "counts_layer": "counts",
            "counts_source": "layers[counts]",
            "counts_validated": False,
            "counts_integer_check": "blockwise",
            "soupx_layer": None,
            "processing_history": [],
            "stage": "02",
        }
        registry[_upstream(root, src, rid).resolve()] = ad
    return {
        "UPSTREAM_RUN_ROOT": str(root), "UPSTREAM_RUNS": runs, "Path": Path, "json": json,
        "pd": pd, "display": lambda _: None, "resume_run": rc.resume_run,
        "sha256_file": rc.sha256_file, "validate_checkpoint": rc.validate_checkpoint,
        "sc": SimpleNamespace(read_h5ad=lambda path: registry[Path(path).resolve()]),
        "sp": sp, "np": np,
    }, registry

def test_json_ast_params_and_prepare_run_placement() -> None:
    for path in [*NB01, NB02]:
        code = [_source(c) for c in _nb(path)["cells"] if c["cell_type"] == "code"]
        ast.parse("\n".join(code), filename=path)
    for path in NB01:
        assert all(name in _params_source(path) for name in ("RUN_ID", "RUN_ROOT", "OUTPUT_FILENAME"))
        final = _cell(path, "# Checkpoint")
        assert "snapshot_effective_parameters(globals()" in final and "collect_runtime_provenance" in final
        earlier = [_source(c) for c in _nb(path)["cells"] if c["cell_type"] == "code" and _source(c) != final]
        assert not any(isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "prepare_run"
                       for source in earlier for n in ast.walk(ast.parse(source)))
    assert all(" in dir()" not in _cell(path, "qc_report = {") for path in (NB01[0], NB01[3]))

@pytest.mark.parametrize("path", NB01)
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

@pytest.mark.parametrize("values", [[], ["Wrong", "Wrong"], ["Kim_2023", "Other"], ["Kim_2023", np.nan]])
def test_stage01_bad_source_or_empty_fails_without_checkpoint(tmp_path: Path, values: list) -> None:
    env = _env01(tmp_path, NB01[0], "failed", n_obs=len(values))
    env["adata"].obs["source_dataset"] = values
    with pytest.raises(RuntimeError, match="Stage 01 FAILED"):
        exec(_cell(NB01[0], "# Checkpoint"), env)
    draft = tmp_path / "runs/failed/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage_status"] == "FAILED" and "checkpoint" not in manifest
    assert not list(draft.glob("*.h5ad"))

def test_stage01_missing_manifest_does_not_claim_run_id(tmp_path: Path) -> None:
    env = _env01(tmp_path, NB01[0], "missing")
    Path(env["MANIFEST_PATH"]).unlink()
    with pytest.raises(FileNotFoundError):
        exec(_cell(NB01[0], "# Checkpoint"), env)
    assert not Path(env["RUN_ROOT"]).exists()

def test_stage01_partial_write_is_removed_and_audited(tmp_path: Path) -> None:
    env = _env01(tmp_path, NB01[0], "write-failed")
    def fail_write(path: Path, **_: object) -> None:
        Path(path).write_bytes(b"partial")
        raise OSError("disk full")
    env["adata"].write_h5ad = fail_write
    with pytest.raises(OSError, match="disk full"):
        exec(_cell(NB01[0], "# Checkpoint"), env)
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


# ---- P0-f: 02 counts 契约 per-source 校验 ----------------------------------------

@pytest.mark.parametrize(("modify", "match"), [
    ("missing_layer", r"layers\['counts'\] 不存在"),
    ("negative", "含负值"),
    ("float_data", "含非整数值"),
    ("shape_mismatch", r"shape.*不一致"),
])
def test_stage02_rejects_invalid_counts_layer(tmp_path: Path, modify: str, match: str) -> None:
    """逐来源 counts 校验：layers['counts'] 缺失/负值/非整数/shape 不对时应抛错."""
    env, registry = _env02(tmp_path)
    for ad in registry.values():
        if modify == "missing_layer":
            ad.layers = {}
        elif modify == "negative":
            ad.layers["counts"] = sp.csr_matrix(np.array([[-1.0, 2.0], [3.0, 4.0]], dtype=np.float32))
        elif modify == "float_data":
            ad.layers["counts"] = sp.csr_matrix(np.array([[1.5, 2.0], [3.7, 4.0]], dtype=np.float32))
        elif modify == "shape_mismatch":
            ad.layers["counts"] = sp.csr_matrix(np.ones((1, 3), dtype=np.float32))
        break
    with pytest.raises(ValueError, match=match):
        exec(_cell(NB02, "# 加载每个 per-dataset h5ad"), env)


def test_stage02_no_preprocessing_done_aggregation() -> None:
    """cell 7207f34e 不含全局 preprocessing_done 聚合逻辑（决策 2：不靠全局状态猜契约）."""
    source = _cell(NB02, "merge_report_v1")
    assert "all_pp_done" not in source, "preprocessing_done 聚合逻辑应已删除"
    assert "preprocessing_done_merged" not in source, "preprocessing_done_merged 字段应已删除"


def test_stage02_checkpoint_imports_validate_expression_contract() -> None:
    """setup cell 导入 validate_expression_contract."""
    setup = _cell(NB02, "# === Setup ===")
    assert "validate_expression_contract" in setup, "setup 应导入 validate_expression_contract"


def test_stage02_checkpoint_writes_expression_contract() -> None:
    """checkpoint cell 写入 stage='02' 的 expression_contract."""
    # 用 'checked_at": "02_merged"' 精确定位 checkpoint cell（setup cell 也有
    # expression_contract 但仅出现在 import 行，不会匹配此 marker）
    checkpoint = _cell(NB02, 'checked_at": "02_merged"')
    assert '"stage": "02"' in checkpoint
    assert '"counts_source": "layers[counts]"' in checkpoint
    assert '"counts_integer_check": "blockwise"' in checkpoint
    assert '"soupx_layer": None' in checkpoint
    assert "validate_expression_contract" in checkpoint
