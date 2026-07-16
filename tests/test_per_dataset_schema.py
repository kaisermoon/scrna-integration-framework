"""test_per_dataset_schema — 测试 per_dataset 输出契约验证函数与常量。

覆盖：valid adata passes、缺失 layer/column/wrong dtype/illegal values/常量拼写。
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from scrna_integration.per_dataset_schema import (
    DOUBLET_CLASS_COL,
    DOUBLET_CONTRACT_KEY,
    DOUBLET_INCLUDE_COL,
    DOUBLET_PREDICTED_COL,
    DOUBLET_SCORE_COL,
    EXPRESSION_CONTRACT_KEY,
    LAYER_COUNTS,
    LAYER_COUNTS_SOUPX,
    QC_N_GENES_COL,
    QC_PCT_MT_COL,
    QC_TOTAL_COUNTS_COL,
    validate_per_dataset_output,
)


def _make_valid_adata(n_obs=20, n_vars=50, seed=42):
    """构建一个完全符合契约的合成 adata，供各测试做基线。"""
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(rng.poisson(5, size=(n_obs, n_vars)).astype(np.float32))
    obs = pd.DataFrame(
        {
            DOUBLET_CLASS_COL: pd.Categorical(
                ["singlet"] * n_obs,
                categories=["singlet", "uncertain", "doublet"],
            ),
            DOUBLET_INCLUDE_COL: [True] * n_obs,
            DOUBLET_PREDICTED_COL: [False] * n_obs,
            DOUBLET_SCORE_COL: np.full(n_obs, np.nan, dtype=np.float64),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"gene_ids": [f"ENSG{str(i).zfill(11)}" for i in range(n_vars)]},
        index=[f"GENE_{i}" for i in range(n_vars)],
    )
    adata = AnnData(X=X, obs=obs, var=var)
    adata.layers[LAYER_COUNTS] = sp.csr_matrix(
        rng.poisson(5, size=(n_obs, n_vars)).astype(np.float32)
    )
    adata.uns[EXPRESSION_CONTRACT_KEY] = {
        "x_scale": "raw_counts",
        "counts_layer": LAYER_COUNTS,
    }
    return adata


# ── 1. valid adata passes ────────────────────────────────────────────────

def test_valid_adata_passes():
    """符合契约的合成 adata 返回 passed=True，无 errors。"""
    adata = _make_valid_adata()
    result = validate_per_dataset_output(adata)
    assert result["passed"] is True
    assert result["errors"] == []


def test_valid_adata_with_doublet_contract():
    """含 doublet_contract 且存在 doublet 细胞时无 errors。"""
    adata = _make_valid_adata(n_obs=10)
    # 注入一个 doublet 细胞
    adata.obs[DOUBLET_CLASS_COL] = pd.Categorical(
        ["singlet"] * 8 + ["uncertain"] + ["doublet"],
        categories=["singlet", "uncertain", "doublet"],
    )
    adata.uns[DOUBLET_CONTRACT_KEY] = {"method": "scrublet", "n_doublet": 1}
    result = validate_per_dataset_output(adata)
    assert result["passed"] is True
    assert result["errors"] == []


# ── 2. missing layers["counts"] ──────────────────────────────────────────

def test_missing_counts_layer():
    """缺少 layers['counts'] → errors 含相应消息。"""
    adata = _make_valid_adata()
    del adata.layers[LAYER_COUNTS]
    result = validate_per_dataset_output(adata)
    assert result["passed"] is False
    assert any(LAYER_COUNTS in e for e in result["errors"])


# ── 3. counts layer 不是稀疏矩阵 ────────────────────────────────────────

def test_counts_layer_not_sparse():
    """layers['counts'] 为 dense → errors。"""
    adata = _make_valid_adata(n_obs=5, n_vars=10)
    adata.layers[LAYER_COUNTS] = adata.layers[LAYER_COUNTS].toarray()
    result = validate_per_dataset_output(adata)
    assert result["passed"] is False
    assert any("非稀疏" in e for e in result["errors"])


# ── 4. counts 含非整数值 ────────────────────────────────────────────────

def test_counts_non_integer():
    """layers['counts'] 含非整数值 → errors。"""
    adata = _make_valid_adata(n_obs=5, n_vars=10)
    c = adata.layers[LAYER_COUNTS].copy()
    c.data = c.data + 0.5  # 破坏整数性
    adata.layers[LAYER_COUNTS] = c
    result = validate_per_dataset_output(adata)
    assert result["passed"] is False
    assert any("非整数" in e for e in result["errors"])


# ── 5. counts 含负值 ────────────────────────────────────────────────────

def test_counts_negative():
    """layers['counts'] 含负值 → errors。"""
    adata = _make_valid_adata(n_obs=5, n_vars=10)
    c = adata.layers[LAYER_COUNTS].copy()
    c.data[0] = -1.0
    adata.layers[LAYER_COUNTS] = c
    result = validate_per_dataset_output(adata)
    assert result["passed"] is False
    assert any("负值" in e for e in result["errors"])


# ── 6. doublet_class 缺失 ────────────────────────────────────────────────

def test_missing_doublet_class_col():
    """obs 中缺少 doublet_class → errors。"""
    adata = _make_valid_adata()
    del adata.obs[DOUBLET_CLASS_COL]
    result = validate_per_dataset_output(adata)
    assert result["passed"] is False
    assert any(DOUBLET_CLASS_COL in e for e in result["errors"])


# ── 7. doublet_class 非法取值 ────────────────────────────────────────────

def test_doublet_class_illegal_value():
    """doublet_class 含非法取值（如 'unknown'）→ errors。"""
    adata = _make_valid_adata(n_obs=5)
    adata.obs[DOUBLET_CLASS_COL] = pd.Categorical(
        ["singlet", "unknown", "doublet", "singlet", "singlet"],
        categories=["singlet", "unknown", "doublet"],
    )
    result = validate_per_dataset_output(adata)
    assert result["passed"] is False
    assert any("非法取值" in e for e in result["errors"])


# ── 8. doublet_class 非 Categorical dtype ────────────────────────────────

def test_doublet_class_not_categorical():
    """doublet_class 为普通 object dtype → errors。"""
    adata = _make_valid_adata(n_obs=5)
    adata.obs[DOUBLET_CLASS_COL] = ["singlet"] * 5
    result = validate_per_dataset_output(adata)
    assert result["passed"] is False
    assert any("Categorical" in e for e in result["errors"])


# ── 9. doublet_include 缺失 ─────────────────────────────────────────────

def test_missing_doublet_include_col():
    """obs 中缺少 doublet_include → errors。"""
    adata = _make_valid_adata()
    del adata.obs[DOUBLET_INCLUDE_COL]
    result = validate_per_dataset_output(adata)
    assert result["passed"] is False
    assert any(DOUBLET_INCLUDE_COL in e for e in result["errors"])


# ── 10. doublet_include 非 bool dtype ────────────────────────────────────

def test_doublet_include_not_bool():
    """doublet_include dtype 为 int 而非 bool → errors。"""
    adata = _make_valid_adata(n_obs=5)
    adata.obs[DOUBLET_INCLUDE_COL] = [1, 1, 0, 1, 1]
    result = validate_per_dataset_output(adata)
    assert result["passed"] is False
    assert any(DOUBLET_INCLUDE_COL in e for e in result["errors"])


# ── 11. expression_contract 缺失 ────────────────────────────────────────

def test_missing_expression_contract():
    """uns 中缺少 expression_contract → errors。"""
    adata = _make_valid_adata()
    del adata.uns[EXPRESSION_CONTRACT_KEY]
    result = validate_per_dataset_output(adata)
    assert result["passed"] is False
    assert any(EXPRESSION_CONTRACT_KEY in e for e in result["errors"])


# ── 12. expression_contract 非 dict ─────────────────────────────────────

def test_expression_contract_not_dict():
    """expression_contract 为字符串而非 dict → errors。"""
    adata = _make_valid_adata()
    adata.uns[EXPRESSION_CONTRACT_KEY] = "not_a_dict"
    result = validate_per_dataset_output(adata)
    assert result["passed"] is False
    assert any("dict" in e for e in result["errors"])


# ── 13. var 缺少基因标识列 ──────────────────────────────────────────────

def test_var_missing_gene_id_cols():
    """var 缺少 gene_ids 和 symbol 列 → warnings。"""
    adata = _make_valid_adata()
    adata.var = pd.DataFrame(index=adata.var_names)
    result = validate_per_dataset_output(adata)
    assert result["passed"] is True  # warning only, not error
    assert any("基因标识" in w for w in result["warnings"])


# ── 14. 常量值拼写 ──────────────────────────────────────────────────────

def test_constant_spellings():
    """验证常量的字符串拼写正确。"""
    assert DOUBLET_CLASS_COL == "doublet_class"
    assert DOUBLET_INCLUDE_COL == "doublet_include"
    assert DOUBLET_PREDICTED_COL == "predicted_doublet"
    assert DOUBLET_SCORE_COL == "doublet_score"
    assert DOUBLET_CONTRACT_KEY == "doublet_contract"
    assert EXPRESSION_CONTRACT_KEY == "expression_contract"
    assert LAYER_COUNTS == "counts"
    assert LAYER_COUNTS_SOUPX == "counts_soupx"


# ── 15. QC 常量 ─────────────────────────────────────────────────────────

def test_qc_constants():
    """验证 QC 指标常量。"""
    assert QC_N_GENES_COL == "n_genes_by_counts"
    assert QC_TOTAL_COUNTS_COL == "total_counts"
    assert QC_PCT_MT_COL == "pct_counts_mt"


# ── 16. 仅 warnings 不阻塞 passed ──────────────────────────────────────

def test_warnings_dont_block_passed():
    """仅有 warnings 没有 errors 时 passed 仍为 True。"""
    adata = _make_valid_adata()
    adata.var = pd.DataFrame(index=adata.var_names)  # 无基因标识列 → warning only
    # 同时确保有 doublet cell 但没有 contract → warning
    adata.obs[DOUBLET_CLASS_COL] = pd.Categorical(
        ["doublet"] + ["singlet"] * (adata.n_obs - 1),
        categories=["singlet", "uncertain", "doublet"],
    )
    result = validate_per_dataset_output(adata)
    assert result["passed"] is True
    assert len(result["warnings"]) >= 1
    assert len(result["errors"]) == 0
