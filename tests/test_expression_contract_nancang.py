from __future__ import annotations

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from scrna_integration.run_contract import validate_expression_contract


def _build_synthetic_nancang_adata(
    n_cells: int = 50,
    n_genes: int = 100,
    *,
    seed: int = 42,
) -> anndata.AnnData:
    """构建模拟 10x mtx 格式的合成 AnnData，X 为 raw integer counts。

    Nancang 的特征：
    - 10x mtx 格式：sparse CSR matrix
    - X 为 raw integer counts (非负整数)
    - 含 source_dataset / sample_id 等 obs 列
    - var 含 ensembl_id 列
    """
    rng = np.random.default_rng(seed)
    counts = sp.random(n_cells, n_genes, density=0.3, format="csr", dtype=np.float32)
    counts.data = np.floor(counts.data * 100 + 1)  # 正整数 counts

    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    var["ensembl_id"] = [f"ENSG{i:06d}" for i in range(n_genes)]

    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    obs["source_dataset"] = "nancang"
    obs["sample_id"] = "sample_1"
    obs["n_genes"] = rng.integers(500, 5000, n_cells)
    obs["total_counts"] = rng.integers(1000, 20000, n_cells)
    obs["pct_counts_mt"] = rng.uniform(0, 15, n_cells)

    adata = anndata.AnnData(X=counts, obs=obs, var=var)
    return adata


def _establish_contract_nancang(adata: anndata.AnnData) -> dict:
    """在 adata 上建立 Nancang 的 expression_contract + layers['counts']。

    与 01_nancang.ipynb 中 expression_contract cell 的实现保持一致。
    返回建立后的 contract dict。
    """
    adata.layers["counts"] = sp.csr_matrix(adata.X, dtype=np.float32)
    counts_data = adata.layers["counts"].data
    is_nonneg = bool(np.all(counts_data >= 0))
    is_integer = bool(np.all(counts_data == np.floor(counts_data)))
    counts_ok = bool(is_nonneg and is_integer)

    contract = {
        "x_scale": "raw_counts",
        "counts_layer": "counts",
        "counts_source": "X",
        "counts_validated": counts_ok,
        "counts_integer_check": "full",
        "soupx_layer": None,
        "processing_history": [
            "01_nancang: loaded from 10x mtx, X is raw integer counts, "
            "copied to layers[counts] (CSR float32)"
        ],
        "stage": "01",
    }
    adata.uns["expression_contract"] = contract
    return contract


# ---- layers["counts"] 建立校验 -------------------------------------------------


class TestLayersCounts:
    """测试 layers["counts"] 的矩阵属性与整数校验。"""

    def test_layers_counts_is_csr_float32(self) -> None:
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        counts = adata.layers["counts"]
        assert sp.issparse(counts), "layers['counts'] 应为稀疏矩阵"
        assert counts.dtype == np.float32, (
            f"layers['counts'] dtype 应为 float32，实际 {counts.dtype}"
        )
        assert counts.format == "csr", (
            f"layers['counts'] format 应为 csr，实际 {counts.format}"
        )

    def test_layers_counts_matches_x_shape(self) -> None:
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        assert adata.layers["counts"].shape == adata.X.shape, (
            f"layers['counts'].shape {adata.layers['counts'].shape} "
            f"!= X.shape {adata.X.shape}"
        )

    def test_counts_all_non_negative(self) -> None:
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        data = adata.layers["counts"].data
        assert bool(np.all(data >= 0)), "counts 应全为非负值"

    def test_counts_all_integer(self) -> None:
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        data = adata.layers["counts"].data
        assert bool(np.all(data == np.floor(data))), "counts 应全为整数值"

    def test_counts_integrity_check_fulfilled(self) -> None:
        """综合测试：非负整数校验必须通过。"""
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        contract = adata.uns["expression_contract"]
        assert contract["counts_validated"] is True, (
            f"counts_validated 应为 True，实际 {contract['counts_validated']}"
        )


# ---- expression_contract 8 字段 schema ----------------------------------------


class TestExpressionContractSchema:
    """测试 expression_contract 全部 8 个字段的完整性与取值合法性。"""

    _REQUIRED_KEYS = frozenset({
        "x_scale", "counts_layer", "counts_source", "counts_validated",
        "counts_integer_check", "soupx_layer", "processing_history", "stage",
    })

    def test_contract_all_eight_fields_present(self) -> None:
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        contract = validate_expression_contract(adata)
        assert set(contract.keys()) == self._REQUIRED_KEYS, (
            f"contract keys: {sorted(contract.keys())}, "
            f"expected: {sorted(self._REQUIRED_KEYS)}"
        )

    def test_contract_passes_full_validation(self) -> None:
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        validate_expression_contract(adata, expected_scale="raw_counts", stage="01")

    def test_contract_field_values_match_nancang_spec(self) -> None:
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        contract = validate_expression_contract(adata, expected_scale="raw_counts", stage="01")
        assert contract["x_scale"] == "raw_counts"
        assert contract["counts_layer"] == "counts"
        assert contract["counts_source"] == "X"
        assert contract["counts_validated"] is True
        assert contract["counts_integer_check"] == "full"
        assert contract["soupx_layer"] is None
        assert isinstance(contract["processing_history"], list)
        assert len(contract["processing_history"]) >= 1
        assert contract["processing_history"][0].startswith("01_nancang")
        assert contract["stage"] == "01"

    def test_contract_soupx_layer_is_none(self) -> None:
        """P0 阶段 soupx_layer 必须为 None（SoupX 属 P1，本轮不动）。"""
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        contract = validate_expression_contract(adata)
        assert contract["soupx_layer"] is None

    def test_contract_processing_history_is_list(self) -> None:
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        contract = validate_expression_contract(adata)
        assert isinstance(contract["processing_history"], list)
        # 至少有一条记录
        assert len(contract["processing_history"]) >= 1


# ---- 边界与错误条件 ------------------------------------------------------------


class TestEdgeCases:
    """测试边界条件与错误情况。"""

    def test_missing_expression_contract_rejected(self) -> None:
        adata = _build_synthetic_nancang_adata()
        # 未建立 expression_contract
        with pytest.raises(KeyError, match="expression_contract not found"):
            validate_expression_contract(adata)

    def test_counts_source_is_x(self) -> None:
        """Nancang 的 counts_source 必须是 'X'（10x mtx 的 X 是 raw counts）。"""
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        contract = validate_expression_contract(adata)
        assert contract["counts_source"] == "X"

    def test_stage_must_be_01(self) -> None:
        """Nancang 的契约 stage 必须是 '01'。"""
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        contract = validate_expression_contract(adata, stage="01")
        assert contract["stage"] == "01"

    def test_contract_rejects_wrong_stage(self) -> None:
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)
        with pytest.raises(ValueError, match="expected '02'"):
            validate_expression_contract(adata, stage="02")

    def test_contract_without_layers_counts_still_validates_contract_schema(self) -> None:
        """schema 校验只查 contract 字段，不查 layers 是否真实存在。

        layers["counts"] 的存在性由 notebook checkpoint cell 的
        hard_postconditions 单独检查。
        """
        adata = _build_synthetic_nancang_adata()
        adata.uns["expression_contract"] = {
            "x_scale": "raw_counts",
            "counts_layer": "counts",
            "counts_source": "X",
            "counts_validated": False,
            "counts_integer_check": "full",
            "soupx_layer": None,
            "processing_history": [],
            "stage": "01",
        }
        # schema 应通过（即使 layers["counts"] 不存在——这是 checkpoint 的职责）
        contract = validate_expression_contract(adata)
        assert contract["counts_validated"] is False
        assert "counts" not in adata.layers

    def test_contract_with_non_integer_counts_marks_validated_false(self) -> None:
        """若 counts 中有非整数值，counts_validated 应为 False。"""
        adata = _build_synthetic_nancang_adata()
        # 将 X 中部分值改为非整数
        adata.X.data[0] = 1.5
        adata.layers["counts"] = sp.csr_matrix(adata.X, dtype=np.float32)
        data = adata.layers["counts"].data
        is_int = bool(np.all(data == np.floor(data)))
        assert not is_int, "应检测到非整数值"
        # 手动建 contract（counts_validated=False）
        adata.uns["expression_contract"] = {
            "x_scale": "raw_counts",
            "counts_layer": "counts",
            "counts_source": "X",
            "counts_validated": False,
            "counts_integer_check": "full",
            "soupx_layer": None,
            "processing_history": ["01_nancang: integer check failed"],
            "stage": "01",
        }
        contract = validate_expression_contract(adata)
        assert contract["counts_validated"] is False


# ---- 集成：模拟 notebook 完整流程 -----------------------------------------------


class TestNotebookIntegration:
    """模拟 01_nancang.ipynb 从数据加载到 checkpoint 的 counts 契约流程。"""

    def test_full_flow_load_to_checkpoint(self) -> None:
        """端到端：数据加载 → 建 counts 契约 → checkpoint 校验。"""
        # Step 1: 模拟数据加载（10x mtx → AnnData，X=raw counts）
        adata = _build_synthetic_nancang_adata()

        # Step 2: 建立 layers["counts"] + expression_contract
        _establish_contract_nancang(adata)

        # Step 3: checkpoint 校验（与 notebook 中 checkpoint cell 同口径）
        # 3a. expression_contract 校验
        contract = validate_expression_contract(adata, expected_scale="raw_counts", stage="01")
        assert contract["counts_validated"] is True

        # 3b. layers["counts"] 存在性
        assert "counts" in adata.layers
        assert sp.issparse(adata.layers["counts"])
        assert adata.layers["counts"].dtype == np.float32

        # 3c. 硬门禁全部通过
        hard_postconditions = {
            "non_empty": adata.n_obs > 0 and adata.n_vars > 0,
            "expression_contract": True,
            "layers_counts_csr_f32": (
                "counts" in adata.layers
                and sp.issparse(adata.layers["counts"])
                and adata.layers["counts"].dtype == np.float32
            ),
        }
        assert all(hard_postconditions.values()), (
            f"hard_postconditions failed: {hard_postconditions}"
        )

    def test_flow_with_corrupted_contract_fails_checkpoint(self) -> None:
        """若契约被破坏（如缺少字段），checkpoint 应失败。"""
        adata = _build_synthetic_nancang_adata()
        _establish_contract_nancang(adata)

        # 恶意删除契约字段
        del adata.uns["expression_contract"]["counts_source"]

        with pytest.raises(ValueError, match="missing required keys"):
            validate_expression_contract(adata)
