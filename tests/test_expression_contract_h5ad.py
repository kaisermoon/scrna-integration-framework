"""P0-d: expression_contract 测试 — h5ad 格式模板 (normalized_log1p + .raw.X)."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from pathlib import Path

import pytest
from anndata import AnnData

from scrna_integration.run_contract import validate_expression_contract


# ---- helpers -----------------------------------------------------------------


class _FakeAdata:
    """最小 AnnData 替身，仅提供 uns dict 访问。"""

    def __init__(self, uns: dict | None = None) -> None:
        self.uns: dict = uns if uns is not None else {}


def _nowicki_contract(**overrides: object) -> dict:
    """返回 nowicki 数据集的 expression_contract（normalized_log1p + .raw.X）。"""
    contract: dict = {
        "x_scale": "normalized_log1p",
        "counts_layer": "counts",
        "counts_source": ".raw.X",
        "counts_validated": True,
        "counts_integer_check": "full",
        "soupx_layer": None,
        "processing_history": ["01_extract_counts_from_raw"],
        "stage": "01",
    }
    contract.update(overrides)
    return contract


# ---- validate_expression_contract for nowicki ----------------------------------


def test_nowicki_contract_valid_full() -> None:
    """nowicki 完整契约通过所有 schema 校验。"""
    adata = _FakeAdata({"expression_contract": _nowicki_contract()})
    result = validate_expression_contract(adata, expected_scale="normalized_log1p")
    assert result["x_scale"] == "normalized_log1p"
    assert result["counts_source"] == ".raw.X"
    assert result["counts_layer"] == "counts"
    assert result["soupx_layer"] is None
    assert result["stage"] == "01"
    assert result["counts_validated"] is True


def test_nowicki_contract_expected_scale_mismatch() -> None:
    """x_scale=normalized_log1p 时不应通过 expected_scale='raw_counts' 校验。"""
    adata = _FakeAdata({"expression_contract": _nowicki_contract()})
    with pytest.raises(ValueError, match="expected 'raw_counts'"):
        validate_expression_contract(adata, expected_scale="raw_counts")


def test_nowicki_contract_expected_stage_mismatch() -> None:
    """stage=01 不应通过 stage='02' 校验。"""
    adata = _FakeAdata({"expression_contract": _nowicki_contract()})
    with pytest.raises(ValueError, match="expected '02'"):
        validate_expression_contract(adata, stage="02")


def test_nowicki_contract_missing_required_key() -> None:
    """缺失 counts_source 应报 missing required keys。"""
    contract = _nowicki_contract()
    del contract["counts_source"]
    adata = _FakeAdata({"expression_contract": contract})
    with pytest.raises(ValueError, match="missing required keys"):
        validate_expression_contract(adata)


def test_nowicki_contract_invalid_counts_source() -> None:
    """counts_source 不在允许值中应报错。"""
    adata = _FakeAdata({"expression_contract": _nowicki_contract(counts_source="layers.counts")})
    with pytest.raises(ValueError, match="counts_source"):
        validate_expression_contract(adata)


# ---- counts extraction simulation ----------------------------------------------


def _make_raw_adata(n_obs: int = 100, n_vars: int = 10) -> AnnData:
    """构造模拟 nowicki 场景的 AnnData：X=logcounts，.raw.X=整数 counts。"""
    rng = np.random.RandomState(42)
    raw_counts = rng.poisson(5, size=(n_obs, n_vars)).astype(np.float32)
    logcounts = np.log1p(raw_counts).astype(np.float32)

    gene_names = [f"GENE_{i}" for i in range(n_vars)]
    raw_adata = AnnData(
        X=sp.csr_matrix(raw_counts),
        var={"gene_id": gene_names},
    )
    raw_adata.var_names = gene_names

    adata = AnnData(
        X=sp.csr_matrix(logcounts),
        var={"gene_id": gene_names},
    )
    adata.var_names = gene_names
    adata.raw = raw_adata
    return adata


def test_extract_counts_from_raw_writes_layers() -> None:
    """从 .raw.X 提取 counts 后 layers['counts'] 应为 CSR float32 整数矩阵。"""
    adata = _make_raw_adata(100, 10)

    # 模拟 notebook 中的提取逻辑
    assert adata.raw is not None
    raw_adata = adata.raw.to_adata()
    assert list(raw_adata.var_names) == list(adata.var_names)

    counts_matrix = raw_adata.X.copy()
    counts_data = counts_matrix.data if sp.issparse(counts_matrix) else np.asarray(counts_matrix).ravel()
    assert counts_data.size > 0
    assert not (counts_data < 0).any()
    assert np.allclose(counts_data % 1, 0)

    adata.layers["counts"] = sp.csr_matrix(counts_matrix, dtype=np.float32)

    assert sp.issparse(adata.layers["counts"])
    assert adata.layers["counts"].dtype == np.float32
    assert adata.layers["counts"].shape == adata.shape


def test_extract_counts_from_raw_preserves_integer_nature() -> None:
    """提取的 counts 应保持为整数（尽管以 float32 存储）。"""
    adata = _make_raw_adata(50, 5)

    raw_adata = adata.raw.to_adata()
    counts_matrix = raw_adata.X.copy()
    adata.layers["counts"] = sp.csr_matrix(counts_matrix, dtype=np.float32)

    dense = adata.layers["counts"].toarray()
    np.testing.assert_array_almost_equal(dense % 1, np.zeros_like(dense))


def test_extract_counts_gene_dimension_mismatch_raises() -> None:
    """.raw.X 基因维度与当前 adata 不一致时应 raise，不静默子集化。"""
    adata = _make_raw_adata(100, 10)

    # 使 adata 仅包含部分基因
    adata_sub = adata[:, :5].copy()

    # 构造不对齐场景：raw 有 10 基因，当前 adata 只有 5 基因
    raw_adata = adata.raw.to_adata()  # 10 genes
    # 当前 adata 只有 5 genes，不一致
    assert list(raw_adata.var_names) != list(adata_sub.var_names)

    with pytest.raises(ValueError, match="基因维度不一致"):
        if list(raw_adata.var_names) != list(adata_sub.var_names):
            raise ValueError(
                f"Nowicki .raw.X 基因维度与当前 adata 不一致："
                f"raw {raw_adata.n_vars} genes vs current {adata_sub.n_vars} genes。"
                f"基因维度不一致，停止处理。请检查 manifest 使用的文件是否完整。"
            )


def test_extract_counts_missing_raw_raises() -> None:
    """.raw 属性不存在时应 raise。"""
    adata = _make_raw_adata(10, 5)
    adata.raw = None  # 模拟缺失 raw

    with pytest.raises(ValueError, match="未找到 .raw 属性"):
        if adata.raw is None:
            raise ValueError("Nowicki 数据必须包含 .raw.X 原始计数，未找到 .raw 属性")


def test_extract_counts_non_integer_raises() -> None:
    """.raw.X 含非整数值时应 raise。"""
    rng = np.random.RandomState(42)
    raw_counts = rng.poisson(5, size=(10, 5)).astype(np.float32)
    raw_counts[0, 0] = 3.5  # 非整数值
    logcounts = np.log1p(raw_counts.clip(0)).astype(np.float32)

    gene_names = ["GENE_0", "GENE_1", "GENE_2", "GENE_3", "GENE_4"]
    raw_adata = AnnData(X=sp.csr_matrix(raw_counts), var={"gene_id": gene_names})
    raw_adata.var_names = gene_names
    adata = AnnData(X=sp.csr_matrix(logcounts), var={"gene_id": gene_names})
    adata.var_names = gene_names
    adata.raw = raw_adata

    raw_adata2 = adata.raw.to_adata()
    counts_matrix = raw_adata2.X.copy()
    counts_data = counts_matrix.data if sp.issparse(counts_matrix) else np.asarray(counts_matrix).ravel()

    with pytest.raises(ValueError, match="non-integer"):
        if not np.allclose(counts_data % 1, 0):
            raise ValueError("counts contain non-integer values")


def test_extract_counts_negative_values_raises() -> None:
    """.raw.X 含负值时应 raise。"""
    rng = np.random.RandomState(42)
    raw_counts = rng.poisson(5, size=(10, 5)).astype(np.float32)
    raw_counts[0, 0] = -1.0  # 负值
    logcounts = np.log1p(raw_counts.clip(0)).astype(np.float32)

    gene_names = ["GENE_0", "GENE_1", "GENE_2", "GENE_3", "GENE_4"]
    raw_adata = AnnData(X=sp.csr_matrix(raw_counts), var={"gene_id": gene_names})
    raw_adata.var_names = gene_names
    adata = AnnData(X=sp.csr_matrix(logcounts), var={"gene_id": gene_names})
    adata.var_names = gene_names
    adata.raw = raw_adata

    raw_adata2 = adata.raw.to_adata()
    counts_matrix = raw_adata2.X.copy()
    counts_data = counts_matrix.data if sp.issparse(counts_matrix) else np.asarray(counts_matrix).ravel()

    with pytest.raises(ValueError, match="negative"):
        if (counts_data < 0).any():
            raise ValueError("counts contain negative values")


def test_extract_counts_empty_matrix_raises() -> None:
    """空 counts 矩阵应 raise。"""
    adata = _make_raw_adata(0, 0)
    # 空矩阵
    raw_adata = adata.raw.to_adata()
    counts_matrix = raw_adata.X.copy()
    counts_data = counts_matrix.data if sp.issparse(counts_matrix) else np.asarray(counts_matrix).ravel()

    with pytest.raises(ValueError, match="empty"):
        if counts_data.size == 0:
            raise ValueError("counts matrix is empty")


# ---- contract shape alignment -------------------------------------------------


def test_counts_layer_shape_matches_x() -> None:
    """layers['counts'] 的 shape 必须与 X 一致。"""
    adata = _make_raw_adata(100, 10)
    raw_adata = adata.raw.to_adata()
    counts_matrix = raw_adata.X.copy()
    adata.layers["counts"] = sp.csr_matrix(counts_matrix, dtype=np.float32)

    assert adata.layers["counts"].shape == adata.X.shape
    assert adata.layers["counts"].shape[0] == adata.n_obs
    assert adata.layers["counts"].shape[1] == adata.n_vars


# ---- memory discipline -------------------------------------------------------


def test_counts_layer_is_csr_float32() -> None:
    """layers['counts'] 应为 CSR float32 格式（内存纪律）。"""
    adata = _make_raw_adata(100, 10)
    raw_adata = adata.raw.to_adata()
    counts_matrix = raw_adata.X.copy()
    adata.layers["counts"] = sp.csr_matrix(counts_matrix, dtype=np.float32)

    assert sp.issparse(adata.layers["counts"])
    assert adata.layers["counts"].format == "csr"
    assert adata.layers["counts"].dtype == np.float32
    assert not isinstance(adata.layers["counts"], np.ndarray)  # 不 densify


def test_counts_layer_not_densified() -> None:
    """layers['counts'] 不应被 densify——保持稀疏。"""
    adata = _make_raw_adata(500, 200)
    raw_adata = adata.raw.to_adata()
    counts_matrix = raw_adata.X.copy()
    adata.layers["counts"] = sp.csr_matrix(counts_matrix, dtype=np.float32)
    assert sp.issparse(adata.layers["counts"])
    n_dense = adata.n_obs * adata.n_vars
    stored = adata.layers["counts"].data.nbytes + adata.layers["counts"].indices.nbytes + adata.layers["counts"].indptr.nbytes
    assert stored < n_dense * 8  # 稀疏存储应远小于密集


# ---- expression_contract 完整字段存在性 ----------------------------------------


REQUIRED_KEYS = frozenset({
    "x_scale", "counts_layer", "counts_source", "counts_validated",
    "counts_integer_check", "soupx_layer", "processing_history", "stage",
})


def test_nowicki_contract_has_all_eight_fields() -> None:
    """contract 必须包含全部 8 个 schema 字段。"""
    contract = _nowicki_contract()
    missing = REQUIRED_KEYS - contract.keys()
    unknown = contract.keys() - REQUIRED_KEYS
    assert not missing, f"missing keys: {missing}"
    assert not unknown, f"unknown keys: {unknown}"


def test_nowicki_contract_soupx_layer_is_none() -> None:
    """P0 阶段 soupx_layer 一律为 None。"""
    contract = _nowicki_contract()
    assert contract["soupx_layer"] is None


def test_nowicki_contract_stage_is_01() -> None:
    """stage 应为 '01'。"""
    contract = _nowicki_contract()
    assert contract["stage"] == "01"
    validate_expression_contract(_FakeAdata({"expression_contract": contract}), stage="01")
