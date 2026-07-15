import ast
import gc
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import scrna_integration.run_contract as rc
ROOT = Path(__file__).parents[1]
NBS = ["notebooks/02_merged.ipynb", "notebooks/03_normalized.ipynb", "notebooks/04_embedded.ipynb"]
def _source(cell: dict) -> str:
    return "".join(value) if isinstance(value := cell.get("source", ""), list) else value
def _cell(path: str, marker: str) -> str:
    return next(_source(cell) for cell in json.loads((ROOT / path).read_text(encoding="utf-8"))["cells"] if marker in _source(cell))
def _params_code_02(path: str) -> str:
    """P1-e 四组化适配：收集 02_merged 的全部 PARAMS code cell（组1-4）。"""
    cells = json.loads((ROOT / path).read_text(encoding="utf-8"))["cells"]
    group_markers = ["组1：数据与版本", "组2：科学参数", "组3：计算参数", "组4：输出与运行标识"]
    parts = []
    for cell in cells:
        if cell["cell_type"] == "code":
            src = _source(cell)
            if any(marker in src for marker in group_markers):
                parts.append(src)
    if len(parts) != 4:
        raise ValueError(f"Expected 4 PARAMS group code cells, found {len(parts)}")
    return "\n".join(parts)
def _cells_all(path: str, marker: str) -> str:
    """P1-e 向后兼容：当参数被拆分到多个 cell 时，拼接所有匹配 cell 的源码。
    用于替代 _cell 的 PARAMS 搜索，使 exec 与断言能跨多 cell 生效。
    """
    nb = json.loads((ROOT / path).read_text(encoding="utf-8"))
    sources = [_source(cell) for cell in nb["cells"] if marker in _source(cell)]
    if not sources:
        raise StopIteration
    return "\n".join(sources)
class _Adata:
    def __init__(self, n_obs: int = 2) -> None:
        self.X = sp.csr_matrix(np.ones((n_obs, 2), dtype=np.float32))
        self.n_obs, self.n_vars = self.X.shape
        self.obs = pd.DataFrame(index=range(n_obs))
        # P1-d：checkpoint cell 会检查 DOUBLET_INCLUDE_KEY 是否在 obs.columns
        self.obs["doublet_include"] = True
        self.var = pd.DataFrame({"highly_variable": [True, False]})
        self.uns: dict = {}
        self.layers = {"counts": sp.csr_matrix(np.ones((n_obs, 2), dtype=np.float32))}

    def write_h5ad(self, path: Path, **_: object) -> None:
        Path(path).write_bytes(b"checkpoint")
def _upstream(root: Path, run_id: str, stage: str) -> Path:
    paths = rc.prepare_run(root, run_id)
    checkpoint = paths.draft_dir / "checkpoint.h5ad"
    checkpoint.write_bytes(b"upstream")
    rc.atomic_write_json(paths.manifest_path, {
        "run_id": run_id, "stage": stage, "stage_status": "SUCCESS",
        "checkpoint": {"path": checkpoint.name, "sha256": rc.sha256_file(checkpoint)},
    })
    return rc.promote_run(paths)
def _base(tmp_path: Path) -> dict:
    names = ("atomic_write_json", "determine_stage_status", "prepare_run", "promote_run",
             "resume_run", "sha256_file", "snapshot_effective_parameters", "validate_checkpoint")
    env = {name: getattr(rc, name) for name in names}
    env.update(Path=Path, json=json, np=np, sp=sp, gc=gc, _root=str(tmp_path),
               collect_runtime_provenance=lambda *_: {"python": "test"})
    return env
def _env02(tmp_path: Path, n_obs: int = 2) -> dict:
    env = _base(tmp_path)
    exec(_params_code_02(NBS[0]), env)  # P1-e 四组化：组1-4 全部 exec
    sources = ["Nancang_2025", "Kim_2023", "Nowicki_2023", "Yue_2024"]
    upstream_inputs = []
    for i, source in enumerate(sources):
        manifest_path, checkpoint_path = tmp_path / f"manifest-{i}.json", tmp_path / f"checkpoint-{i}.h5ad"
        manifest_path.write_text("{}", encoding="utf-8"); checkpoint_path.write_bytes(b"upstream")
        upstream_inputs.append({"source_dataset": source, "run_id": f"run-{i}", "stage": "01_qcd", "manifest_path": str(manifest_path),
            "manifest_sha256": rc.sha256_file(manifest_path), "checkpoint_path": str(checkpoint_path),
            "checkpoint_sha256": rc.sha256_file(checkpoint_path)})
    env["validate_expression_contract"] = rc.validate_expression_contract
    # P1-d：checkpoint cell 引用的 doublet 集成变量
    env.update(adata=_Adata(n_obs), upstream_inputs=upstream_inputs, PER_DATASET_PATHS=sources, HOUSEKEEPING_GENES=["ACTB", "GAPDH"],
               RUN_ID="stage02", RUN_ROOT=str(tmp_path / "runs"),
               OUTPUT_FILENAME="02.h5ad",
               doublet_inclusion_report={"n_before": n_obs, "n_after": n_obs, "n_excluded": 0},
               doublet_needs_review=False, flagged_samples=[])
    return env
def _env03(tmp_path: Path, stage: str = "02_merged") -> dict:
    root = tmp_path / "upstream"
    checkpoint = _upstream(root, "stage02", stage)
    adata = _Adata()
    adata.uns["expression_contract"] = {
        "x_scale": "raw_counts", "counts_layer": "counts",
        "counts_source": "layers[counts]", "counts_validated": True,
        "counts_integer_check": "blockwise", "soupx_layer": None,
        "processing_history": [], "stage": "02",
    }
    env = _base(tmp_path)
    exec(_cells_all(NBS[1], "# === PARAMS ==="), env)
    env.update(UPSTREAM_RUN_ROOT=str(root), UPSTREAM_RUN_ID="stage02", RUN_ID="stage03",
               RUN_ROOT=str(tmp_path / "runs"), OUTPUT_FILENAME="03.h5ad",
               sc=SimpleNamespace(read_h5ad=lambda path: adata if Path(path) == checkpoint else None),
               excluded_counts={}, _n_genes_values=[3000], _flavor_values=["seurat"],
               validate_expression_contract=rc.validate_expression_contract)
    return env
def test_json_ast_params_and_prepare_run_is_deferred() -> None:
    for path in NBS:
        notebook = json.loads((ROOT / path).read_text(encoding="utf-8"))
        code = [_source(cell) for cell in notebook["cells"] if cell["cell_type"] == "code"]
        ast.parse("\n".join(code), filename=path)
        if path != NBS[2]:
            # P1-e：02_merged(NBS[0]) 四组化用 _params_code_02；03(NBS[1]) 四组化用 _cells_all；其余单 cell
            if path == NBS[0]:
                params = _params_code_02(path)
            elif path == NBS[1]:
                params = _cells_all(path, "# === PARAMS ===")
            else:
                params = _cell(path, "# === PARAMS ===")
            assert all(name in params for name in ("RUN_ID", "RUN_ROOT", "OUTPUT_FILENAME"))
            final = _cell(path, "# Checkpoint" if path == NBS[0] else "# === 参数记录 + Checkpoint ===")
            earlier = [source for source in code if source != final]
            assert not any(isinstance(node, ast.Call) and getattr(node.func, "id", None) == "prepare_run"
                           for source in earlier for node in ast.walk(ast.parse(source)))
@pytest.mark.parametrize("case", ["success", "gate", "provenance", "write"])
def test_stage02_promotion_and_failure_audit(tmp_path: Path, case: str) -> None:
    env = _env02(tmp_path, 0 if case == "gate" else 2)
    adata_ref = env["adata"]
    if case == "provenance":
        env["upstream_inputs"][0].pop("checkpoint_sha256")
    if case == "write":
        def fail(path: Path, **_: object) -> None:
            Path(path).write_bytes(b"partial")
            raise OSError("disk full")
        env["adata"].write_h5ad = fail
    error = pytest.raises(RuntimeError if case in {"gate", "provenance"} else OSError) if case != "success" else None
    if error:
        with error:
            exec(_cell(NBS[0], "# Checkpoint"), env)
        manifest_path = tmp_path / "runs/stage02/draft/manifest.json"
        manifest = json.loads(manifest_path.read_text())
        assert manifest["stage_status"] == "FAILED" and "checkpoint" not in manifest
        assert not list(manifest_path.parent.glob("*.h5ad"))
    else:
        exec(_cell(NBS[0], "# Checkpoint"), env)
        manifest_path = tmp_path / "runs/stage02/promoted/manifest.json"
        manifest = json.loads(manifest_path.read_text())
        assert manifest["stage"] == "02_merged" and len(manifest["inputs"]) == 4
        assert {"BATCH_KEY", "HOUSEKEEPING_GENES"} <= manifest["effective_parameters"].keys() and "PER_DATASET_PATHS" not in manifest["effective_parameters"]
        assert adata_ref.uns["upstream_inputs"] == env["upstream_inputs"]
        assert "adata" not in env
        assert rc.validate_checkpoint(manifest_path).is_file()
@pytest.mark.parametrize("case", ["success", "write"])
def test_stage03_verified_input_and_promotion(tmp_path: Path, case: str) -> None:
    env = _env03(tmp_path)
    exec(_cell(NBS[1], "# === 加载上游 ==="), env)
    # P0-g: 为 normalize cell 提供 mock sc.pp（免真实 scanpy），
    # 补 exec counts-layer 与 normalize-code cells，
    # 使 _upstream_contract 和 _counts_integrity_checked 在 checkpoint cell 前已定义
    env["sc"].pp = SimpleNamespace(
        normalize_total=lambda adata, target_sum: None,
        log1p=lambda adata: None,
    )
    exec(_cell(NBS[1], "从上游验证 layers"), env)
    exec(_cell(NBS[1], "# === 标准化分支 ==="), env)
    adata_ref = env["adata"]
    if case == "write":
        def fail(path: Path, **_: object) -> None:
            Path(path).write_bytes(b"partial")
            raise OSError("disk full")
        env["adata"].write_h5ad = fail
        with pytest.raises(OSError, match="disk full"):
            exec(_cell(NBS[1], "# === 参数记录 + Checkpoint ==="), env)
        draft = tmp_path / "runs/stage03/draft"
        manifest = json.loads((draft / "manifest.json").read_text())
        assert manifest["stage_status"] == "FAILED" and not list(draft.glob("*.h5ad"))
        return
    exec(_cell(NBS[1], "# === 参数记录 + Checkpoint ==="), env)
    manifest_path = tmp_path / "runs/stage03/promoted/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["stage"] == "03_normalized" and manifest["inputs"][0]["stage"] == "02_merged"
    assert "NORMALIZATION_METHOD" in manifest["effective_parameters"] and "UPSTREAM_CHECKPOINT" not in manifest["effective_parameters"]
    assert adata_ref.uns["upstream_inputs"] == manifest["inputs"]
    assert "adata" not in env
@pytest.mark.parametrize(("case", "match"), [("stage", "stage"), ("hash", "hash")])
def test_stage03_rejects_invalid_upstream(tmp_path: Path, case: str, match: str) -> None:
    env = _env03(tmp_path, "wrong" if case == "stage" else "02_merged")
    if case == "hash":
        (Path(env["UPSTREAM_RUN_ROOT"]) / "stage02/promoted/checkpoint.h5ad").write_bytes(b"bad")
    with pytest.raises(ValueError, match=match):
        exec(_cell(NBS[1], "# === 加载上游 ==="), env)
@pytest.mark.parametrize(("case", "match"), [
    ("ok", None), ("warnings", None), ("stage", "stage"), ("hash", "hash"),
    ("missing_status", "stage_status"), ("failed", "stage_status"), ("needs_review", "stage_status"),
    ("unaccepted_warnings", "warning_acceptance"),
])
def test_stage04_reads_only_verified_stage03(tmp_path: Path, case: str, match: str | None) -> None:
    root = tmp_path / "runs"
    checkpoint = _upstream(root, "stage03", "wrong" if case == "stage" else "03_normalized")
    if case == "hash":
        checkpoint.write_bytes(b"bad")
    manifest_path = root / "stage03/promoted/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if case == "missing_status":
        manifest.pop("stage_status")
    elif case in {"failed", "needs_review"}:
        manifest["stage_status"] = case.upper()
    elif case in {"warnings", "unaccepted_warnings"}:
        manifest["stage_status"] = "SUCCESS_WITH_WARNINGS"
        manifest["warning_acceptance"] = ({"accepted_by": "researcher", "accepted_at": "2026-07-13T12:00:00Z"} if case == "warnings" else {})
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    adata = _Adata()
    env = {"UPSTREAM_RUN_ROOT": str(root), "UPSTREAM_RUN_ID": "stage03", "Path": Path,
           "json": json, "resume_run": rc.resume_run, "sha256_file": rc.sha256_file,
           "validate_checkpoint": rc.validate_checkpoint, "sp": sp,
           "sc": SimpleNamespace(read_h5ad=lambda path: adata if Path(path) == checkpoint else None)}
    code = (source := _cell(NBS[2], "# === 只读取 promoted Stage 03 上游 ==="))[source.index("# === 只读取 promoted Stage 03 上游 ==="):source.index("# === promoted Stage 03 上游读取完成 ===")]
    if match:
        with pytest.raises(ValueError, match=match):
            exec(code, env)
    else:
        exec(code, env)
        assert env["adata"] is adata and env["upstream_input"]["stage"] == "03_normalized"
