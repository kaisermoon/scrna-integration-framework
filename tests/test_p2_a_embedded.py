"""P2-a-embedded（UX-1 + 决策4）叶子测试：04_embedded.ipynb 的方法勾选、
scVI 严格六项校验、per-method 状态与研究者显式 selected_embedding 决策 cell。

测试策略（遵循 notebook 测试纪律）：直接 exec cell 源码，用轻量 mock adata 验证
契约逻辑与决策门禁，不跑真实 scanpy/scVI 重计算。
"""
import ast
import json
import re
from enum import Enum
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

ROOT = Path(__file__).parents[1]
NOTEBOOK = ROOT / "notebooks/04_embedded.ipynb"


class MethodStatus(str, Enum):
    SUCCESS = "success"
    SKIPPED_BY_USER = "skipped_by_user"
    UNAVAILABLE = "unavailable"
    FAILED = "failed"


def _source(cell: dict) -> str:
    value = cell.get("source", "")
    return "".join(value) if isinstance(value, list) else value


def _find_cell(marker: str) -> str:
    """在 notebook 的 code cell 中按 marker 子串查找，返回源码。"""
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code" and marker in _source(cell):
            return _source(cell)
    raise LookupError(f"marker not found: {marker!r}")


def _base_env() -> dict:
    """exec PARAMS 拿到全部默认参数，注入 MethodStatus 等运行时依赖。"""
    env = {
        "pd": pd, "np": np, "sp": sp, "re": re, "MethodStatus": MethodStatus,
        "UPSTREAM_RUN_ROOT": "x", "UPSTREAM_RUN_ID": "x",
    }
    exec(_find_cell("# === PARAMS ==="), env)
    return env


def _simple_adata(n_obs=4, n_vars=10, **obs_cols) -> SimpleNamespace:
    """构建轻量 mock adata（SimpleNamespace），不拉真实 anndata 开销。"""
    obs = pd.DataFrame(obs_cols, index=range(n_obs))
    return SimpleNamespace(
        obs=obs, uns={}, obsm={}, varm={}, layers={},
        n_obs=n_obs, n_vars=n_vars,
        var_names=pd.Index([f"gene_{i}" for i in range(n_vars)]),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)]),
        X=sp.csr_matrix(np.ones((n_obs, n_vars), dtype=np.float32)),
    )


# ============================================================================
# 静态结构：AST + marker 保留 + 状态标记纪律（P0-i）
# ============================================================================

def test_all_code_cells_parse_and_markers_preserved() -> None:
    """全 code cell ast.parse 合法 + UX-1/决策4 关键 token 在位。"""
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    code = [_source(c) for c in notebook["cells"] if c["cell_type"] == "code"]
    joined = "\n".join(code)
    ast.parse(joined, filename=str(NOTEBOOK))

    for token in (
        # UX-1 per-method flags
        "PCA_ENABLED", "HARMONY_ENABLED", "SCVI_ENABLED",
        "SCANVI_ENABLED", "CENSUS_ENABLED", "SCCRAFT_ENABLED",
        # UX-1 decision cell
        "SELECTED_EMBEDDING", "SELECTION_RATIONALE",
        # Decision 4
        "_counts_valid", "_counts_failure_reasons",
        "_method_status",
    ):
        assert token in joined, token


def test_new_cells_use_explicit_completed_markers() -> None:
    """P0-i：新增 cell 用显式 _xxx_completed，不用 dir() 判存在性。"""
    # 决策4 scVI 严格校验 cell 使用显式 _counts_valid 标记
    valid = _find_cell("# === scVI/scANVI 输入严格校验（决策 4）===")
    assert "_counts_valid = False" in valid
    assert "_counts_valid = all(_all_checks)" in valid
    assert "in dir()" not in valid

    # selected_embedding 决策 cell
    decision = _find_cell("# === 研究者显式选择 selected_embedding")
    assert "_decision_completed = False" in decision
    assert "_decision_completed = True" in decision
    assert "in dir()" not in decision


def test_params_defines_ux1_and_decision_defaults() -> None:
    """PARAMS 中 UX-1 per-method flags 与决策字段默认值正确。"""
    env = _base_env()

    # Per-method boolean flags
    assert env["PCA_ENABLED"] is True
    assert env["HARMONY_ENABLED"] is True
    assert env["SCVI_ENABLED"] is True
    assert env["SCANVI_ENABLED"] is True
    assert env["CENSUS_ENABLED"] is False
    assert env["SCCRAFT_ENABLED"] is False

    # Decision cell defaults
    assert env["SELECTED_EMBEDDING"] is None
    assert env["SELECTION_RATIONALE"] == ""


# ============================================================================
# 决策4：scVI 严格六项校验
# ============================================================================

def test_scvi_validation_all_six_checks_present() -> None:
    """决策4 六项校验全部在代码中出现。"""
    source = _find_cell("# === scVI/scANVI 输入严格校验（决策 4）===")
    checks = [
        "shape",           # 1. shape 一致
        "gene",            # 2. 基因顺序
        "finite",          # 3. 有限非负
        "integer",         # 4. 近整数
        "library",         # 5. 每批次有效文库大小
        "contract",        # 6. 契约元数据
    ]
    for chk in checks:
        assert chk in source.lower(), f"Missing check {chk!r}"


def test_scvi_validation_block_on_failure() -> None:
    """决策4 任一校验失败必须 raise 阻断，绝不静默降级 PCA。"""
    all_code = [_source(c) for c in
                json.loads(NOTEBOOK.read_text(encoding="utf-8"))["cells"]
                if c["cell_type"] == "code"]
    joined = "\n".join(all_code)

    # 严格校验失败应有阻断逻辑
    assert "raise" in joined

    # scVI cell 中 _counts_valid=False 导致 UNAVAILABLE
    scvi_src = _find_cell("# scVI：setup_anndata + 训练 + 提取潜表示。")
    assert "MethodStatus.UNAVAILABLE" in scvi_src

    # scANVI cell 同样依赖 _counts_valid
    scanvi_src = _find_cell("# scANVI：scVI 的半监督变体，利用细胞类型标签。")
    assert "MethodStatus.UNAVAILABLE" in scanvi_src


# ============================================================================
# UX-1：方法勾选 + 前置条件检查
# ============================================================================

def test_method_eligibility_preflight_cell_exists() -> None:
    """UX-1 preflight cell 存在并列出所有方法状态。"""
    source = _find_cell("# === 方法可运行性预检（UX-1）===")
    for method in ["pca", "harmony", "scvi", "scanvi", "census", "sccraft"]:
        assert method in source.lower(), f"Preflight missing: {method}"


def test_pca_cell_uses_pca_enabled() -> None:
    """PCA cell 使用 PCA_ENABLED 而非旧列表。"""
    source = _find_cell("# 在高可变基因（HVG）上做主成分分析。")
    assert "PCA_ENABLED" in source
    assert '"pca" in EMBEDDING_METHODS' not in source


def test_harmony_cell_uses_harmony_enabled() -> None:
    """Harmony cell 使用 HARMONY_ENABLED 而非旧列表。"""
    source = _find_cell("# Harmony 在 PCA 嵌入上做批次整合。")
    assert "HARMONY_ENABLED" in source
    assert '"harmony" in EMBEDDING_METHODS' not in source


def test_scvi_cell_uses_scvi_enabled() -> None:
    """scVI cell 使用 SCVI_ENABLED 而非旧列表。"""
    source = _find_cell("# scVI：setup_anndata + 训练 + 提取潜表示。")
    assert "SCVI_ENABLED" in source
    assert '"scvi" in EMBEDDING_METHODS' not in source


def test_scanvi_cell_uses_scanvi_enabled() -> None:
    """scANVI cell 使用 SCANVI_ENABLED 而非旧列表。"""
    source = _find_cell("# scANVI：scVI 的半监督变体，利用细胞类型标签。")
    assert "SCANVI_ENABLED" in source
    assert '"scanvi" in EMBEDDING_METHODS' not in source


def test_sccraft_cell_uses_sccraft_enabled() -> None:
    """scCRAFT cell 使用 SCCRAFT_ENABLED 而非旧列表。"""
    source = _find_cell("# scCRAFT：在独立 counts 副本上训练，结果对齐写回主 adata。")
    assert "SCCRAFT_ENABLED" in source
    assert '"sccraft" in EMBEDDING_METHODS' not in source


def test_census_cell_mentions_census_enabled() -> None:
    """Census cell 提及 CENSUS_ENABLED。"""
    source = _find_cell("# # === cellxgene_census 预训练 scVI（已注释） ===")
    assert "CENSUS_ENABLED" in source


# ============================================================================
# UX-1：selected_embedding 决策 cell
# ============================================================================

def _decision_env(selected=None, embedding_keys=None) -> dict:
    """构建 selected_embedding 决策 cell 的执行环境。"""
    env = _base_env()
    env["SELECTED_EMBEDDING"] = selected
    env["SELECTION_RATIONALE"] = "scVI 批次校正效果最优"
    if embedding_keys is None:
        embedding_keys = ["X_pca", "X_pca_harmony", "X_scVI"]
    env["adata"] = _simple_adata()
    env["adata"].obsm = {k: np.random.randn(4, 10) for k in embedding_keys}
    env["_method_status"] = {
        "pca": MethodStatus.SUCCESS,
        "harmony": MethodStatus.SUCCESS,
        "scvi": MethodStatus.SUCCESS,
    }
    return env


def test_decision_none_keeps_needs_review_no_raise() -> None:
    """未选时不 raise，不自动赋值——护 Run-All / CI。"""
    env = _decision_env(None)
    exec(_find_cell("# === 研究者显式选择 selected_embedding"), env)
    assert "selected_embedding" not in env["adata"].uns


def test_decision_valid_writes_uns() -> None:
    """选择有效嵌入时写入 adata.uns。"""
    env = _decision_env("X_scVI")
    exec(_find_cell("# === 研究者显式选择 selected_embedding"), env)
    assert env["adata"].uns["selected_embedding"] == "X_scVI"
    assert env["adata"].uns["selection_rationale"] == "scVI 批次校正效果最优"


def test_decision_key_not_in_obsm_raises() -> None:
    """选择的 key 不在 obsm 中必须 raise ValueError。"""
    env = _decision_env("X_nonexistent")
    with pytest.raises(ValueError, match=r"不在 adata\.obsm"):
        exec(_find_cell("# === 研究者显式选择 selected_embedding"), env)


def test_decision_does_not_auto_assign_selected_key() -> None:
    """决策 cell 绝不根据指标自动给 SELECTED_EMBEDDING 赋计算值。"""
    source = _find_cell("# === 研究者显式选择 selected_embedding")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for tgt in node.targets:
                assert not (isinstance(tgt, ast.Name)
                            and tgt.id == "SELECTED_EMBEDDING"), \
                    "决策 cell 不得自动赋值 SELECTED_EMBEDDING"


def test_selected_embedding_not_auto_filled_in_any_cell() -> None:
    """全 notebook 除 PARAMS 外任何 cell 不得自动给 SELECTED_EMBEDDING 赋值。"""
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    all_code = [_source(c) for c in notebook["cells"]
                if c["cell_type"] == "code"]
    for cell_source in all_code:
        # PARAMS cell 中 SELECTED_EMBEDDING = None 是用户默认值，不等于自动赋值
        if "# === PARAMS ===" in cell_source:
            continue
        tree = ast.parse(cell_source)
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                for tgt in node.targets:
                    assert not (isinstance(tgt, ast.Name)
                                and tgt.id == "SELECTED_EMBEDDING"), \
                        "SELECTED_EMBEDDING 不得被自动赋值"


# ============================================================================
# 状态机：per-method status + SUCCESS_WITH_WARNINGS
# ============================================================================

def test_per_method_status_initialized() -> None:
    """_method_status 在 PARAMS 中预初始化为空 dict。"""
    source = _find_cell("# === PARAMS ===")
    assert "_method_status" in source


def test_success_with_warnings_in_checkpoint() -> None:
    """checkpoint cell 支持 SUCCESS_WITH_WARNINGS 路由。"""
    source = _find_cell("# === Stage 04 draft checkpoint")
    assert "SUCCESS_WITH_WARNINGS" in source


def test_checkpoint_includes_selected_embedding() -> None:
    """manifest 包含 selected_embedding 和 selection_rationale。"""
    source = _find_cell("# === Stage 04 draft checkpoint")
    assert "selected_embedding" in source
    assert "selection_rationale" in source


def test_checkpoint_includes_counts_validation() -> None:
    """manifest 包含决策4 counts 校验结果。"""
    source = _find_cell("# === Stage 04 draft checkpoint")
    assert "_counts_valid" in source
    assert "_counts_failure_reasons" in source


# ============================================================================
# redline 7 守护：新增 uns 列式结构必须能安全 h5ad round-trip
# ============================================================================

def test_new_uns_keys_survive_h5ad_roundtrip(tmp_path: Path) -> None:
    """P2-a 新增 uns 键在 h5ad round-trip 中不崩溃。"""
    adata_obj = ad.AnnData(sp.csr_matrix(np.ones((4, 2), dtype=np.float32)))
    adata_obj.uns["selected_embedding"] = "X_scVI"
    adata_obj.uns["selection_rationale"] = "批次混合与生物学分离兼优"
    adata_obj.uns["scvi_validation"] = {
        "counts_valid": True,
        "failure_reasons": [],
        "checks": {
            "shape": True, "gene_order": True, "finite_nonnegative": True,
            "near_integer": True, "library_size": True, "contract_metadata": True,
        },
    }
    path = tmp_path / "roundtrip.h5ad"
    adata_obj.write_h5ad(path)
    restored = ad.read_h5ad(path)
    assert restored.uns["selected_embedding"] == "X_scVI"
    assert restored.uns["selection_rationale"] == "批次混合与生物学分离兼优"
    assert restored.uns["scvi_validation"]["checks"]["shape"] is True
