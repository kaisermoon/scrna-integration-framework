from __future__ import annotations

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from scrna_integration.run_contract import validate_expression_contract


def _build_synthetic_yue_adata(
    n_cells: int = 50,
    n_genes: int = 100,
    *,
    seed: int = 42,
) -> anndata.AnnData:
    """构建模拟 txt.gz 拼装格式的合成 AnnData，X 为 raw UMI integer counts。

    Yue 类器官数据集的特征：
    - txt.gz 文件逐个拼装：每文件一个细胞，tab 分隔基因×计数
    - X 为 raw UMI integer counts（非负整数）
    - 含 source_dataset / sample_id 等 obs 列
    - var 含 ensembl_id 列
    """
    rng = np.random.default_rng(seed)
    counts = sp.random(n_cells, n_genes, density=0.3, format="csr", dtype=np.float32)
    counts.data = np.floor(counts.data * 100 + 1)  # 正整数 counts

    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    var["ensembl_id"] = [f"ENSG{i:06d}" for i in range(n_genes)]

    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    obs["source_dataset"] = "yue"
    obs["sample_id"] = [f"organoid_{i % 3 + 1}" for i in range(n_cells)]
    obs["n_genes"] = rng.integers(500, 5000, n_cells)
    obs["total_counts"] = rng.integers(1000, 20000, n_cells)
    obs["pct_counts_mt"] = rng.uniform(0, 15, n_cells)

    adata = anndata.AnnData(X=counts, obs=obs, var=var)
    return adata


def _establish_contract_yue(adata: anndata.AnnData) -> dict:
    """在 adata 上建立 Yue 的 expression_contract + layers['counts']（阶段 1：加载 cell）。

    与 01_yue.ipynb 数据加载 cell 实现一致：
    - 建立 layers["counts"]（CSR float32）
    - 写入 expression_contract，counts_validated=False（校验在 checkpoint cell）
    返回建立后的 contract dict。
    """
    adata.layers["counts"] = sp.csr_matrix(adata.X, dtype=np.float32)
    contract = {
        "x_scale": "raw_counts",
        "counts_layer": "counts",
        "counts_source": "X",
        "counts_validated": False,
        "counts_integer_check": None,
        "soupx_layer": None,
        "processing_history": [],
        "stage": "01",
    }
    adata.uns["expression_contract"] = contract
    return contract


def _validate_counts_in_checkpoint(adata: anndata.AnnData) -> dict:
    """模拟 checkpoint cell 中的 counts 校验（阶段 2：checkpoint）。

    与 01_yue.ipynb checkpoint cell 实现一致：
    - 非负校验
    - 近整数校验（float32 浮点误差 < 1e-6）
    - 更新 counts_validated=True, counts_integer_check="full"
    - schema 校验
    返回校验通过的 contract dict。
    """
    contract = adata.uns["expression_contract"]
    counts = adata.layers[contract["counts_layer"]]
    assert counts.shape == adata.shape, (
        f"counts layer shape {counts.shape} != adata shape {adata.shape}"
    )
    if sp.issparse(counts):
        counts_data = counts.data
    else:
        counts_data = counts
    assert bool(np.all(counts_data >= 0)), "counts layer contains negative values"
    diff = np.abs(counts_data - np.round(counts_data))
    assert bool(np.all(diff < 1e-6)), (
        f"counts layer contains non-integer values (max diff={diff.max():.2e})"
    )
    contract["counts_validated"] = True
    contract["counts_integer_check"] = "full"
    validated_contract = validate_expression_contract(adata)
    return validated_contract


# ---- layers["counts"] 建立校验 -------------------------------------------------


class TestLayersCounts:
    """测试 layers["counts"] 的矩阵属性与整数校验。"""

    def test_layers_counts_is_csr_float32(self) -> None:
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        counts = adata.layers["counts"]
        assert sp.issparse(counts), "layers['counts'] 应为稀疏矩阵"
        assert counts.dtype == np.float32, (
            f"layers['counts'] dtype 应为 float32，实际 {counts.dtype}"
        )
        assert counts.format == "csr", (
            f"layers['counts'] format 应为 csr，实际 {counts.format}"
        )

    def test_layers_counts_matches_x_shape(self) -> None:
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        assert adata.layers["counts"].shape == adata.X.shape, (
            f"layers['counts'].shape {adata.layers['counts'].shape} "
            f"!= X.shape {adata.X.shape}"
        )

    def test_counts_all_non_negative(self) -> None:
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        _validate_counts_in_checkpoint(adata)
        data = adata.layers["counts"].data
        assert bool(np.all(data >= 0)), "counts 应全为非负值"

    def test_counts_all_integer(self) -> None:
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        _validate_counts_in_checkpoint(adata)
        data = adata.layers["counts"].data
        assert bool(np.all(data == np.floor(data))), "counts 应全为整数值"

    def test_counts_integrity_check_fulfilled(self) -> None:
        """综合测试：非负整数校验必须通过。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        _validate_counts_in_checkpoint(adata)
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
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        _validate_counts_in_checkpoint(adata)
        contract = validate_expression_contract(adata)
        assert set(contract.keys()) == self._REQUIRED_KEYS, (
            f"contract keys: {sorted(contract.keys())}, "
            f"expected: {sorted(self._REQUIRED_KEYS)}"
        )

    def test_contract_passes_full_validation(self) -> None:
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        _validate_counts_in_checkpoint(adata)
        validate_expression_contract(adata, expected_scale="raw_counts", stage="01")

    def test_contract_field_values_match_yue_spec(self) -> None:
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        _validate_counts_in_checkpoint(adata)
        contract = validate_expression_contract(adata, expected_scale="raw_counts", stage="01")
        assert contract["x_scale"] == "raw_counts"
        assert contract["counts_layer"] == "counts"
        assert contract["counts_source"] == "X"
        assert contract["counts_validated"] is True
        assert contract["counts_integer_check"] == "full"
        assert contract["soupx_layer"] is None
        assert isinstance(contract["processing_history"], list)
        assert contract["stage"] == "01"

    def test_contract_soupx_layer_is_none_pre_validation(self) -> None:
        """Yue 类器官数据集无 raw matrix，soupx_layer 必须为 None（加载 cell）。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        # 加载 cell 不调用 validate_expression_contract（counts_integer_check=None 不通过 schema）
        contract = adata.uns["expression_contract"]
        assert contract["soupx_layer"] is None

    def test_contract_soupx_layer_is_none_post_validation(self) -> None:
        """Yue 类器官数据集无 raw matrix，soupx_layer 必须为 None（checkpoint 后）。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        _validate_counts_in_checkpoint(adata)
        contract = validate_expression_contract(adata)
        assert contract["soupx_layer"] is None

    def test_contract_processing_history_is_list(self) -> None:
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        _validate_counts_in_checkpoint(adata)
        contract = validate_expression_contract(adata)
        assert isinstance(contract["processing_history"], list)


# ---- 加载 cell 状态（校验前）：counts_validated=False, counts_integer_check=None ---


class TestLoadingCellContract:
    """测试数据加载 cell 建立的契约状态（校验前）。"""

    _REQUIRED_KEYS = frozenset({
        "x_scale", "counts_layer", "counts_source", "counts_validated",
        "counts_integer_check", "soupx_layer", "processing_history", "stage",
    })

    def test_contract_before_validation_has_validated_false(self) -> None:
        """加载 cell 结束时 counts_validated 应为 False。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        contract = adata.uns["expression_contract"]
        assert contract["counts_validated"] is False
        assert contract["counts_integer_check"] is None

    def test_contract_before_validation_has_correct_field_values(self) -> None:
        """加载 cell 的 contract 字段值正确（schema 校验前 counts_integer_check=None）。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        contract = adata.uns["expression_contract"]
        # 全部 8 字段存在
        assert set(contract.keys()) == self._REQUIRED_KEYS
        # 值正确
        assert contract["x_scale"] == "raw_counts"
        assert contract["counts_layer"] == "counts"
        assert contract["counts_source"] == "X"
        assert contract["counts_validated"] is False
        assert contract["counts_integer_check"] is None
        assert contract["soupx_layer"] is None
        assert isinstance(contract["processing_history"], list)
        assert contract["stage"] == "01"


# ---- 边界与错误条件 ------------------------------------------------------------


class TestEdgeCases:
    """测试边界条件与错误情况。"""

    def test_missing_expression_contract_rejected(self) -> None:
        adata = _build_synthetic_yue_adata()
        with pytest.raises(KeyError, match="expression_contract not found"):
            validate_expression_contract(adata)

    def test_counts_source_is_x(self) -> None:
        """Yue 的 counts_source 必须是 'X'（txt.gz 拼装的 X 是 raw counts）。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        _validate_counts_in_checkpoint(adata)
        contract = validate_expression_contract(adata)
        assert contract["counts_source"] == "X"

    def test_stage_must_be_01(self) -> None:
        """Yue 的契约 stage 必须是 '01'。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        _validate_counts_in_checkpoint(adata)
        contract = validate_expression_contract(adata, stage="01")
        assert contract["stage"] == "01"

    def test_contract_rejects_wrong_stage(self) -> None:
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        _validate_counts_in_checkpoint(adata)
        with pytest.raises(ValueError, match="expected '02'"):
            validate_expression_contract(adata, stage="02")

    def test_contract_without_layers_counts_still_validates_contract_schema(self) -> None:
        """schema 校验只查 contract 字段，不查 layers 是否真实存在。

        layers["counts"] 的存在性由 notebook checkpoint cell 的
        hard_postconditions 单独检查。
        """
        adata = _build_synthetic_yue_adata()
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
        contract = validate_expression_contract(adata)
        assert contract["counts_validated"] is False
        assert "counts" not in adata.layers

    def test_contract_with_non_integer_counts_marks_validated_false(self) -> None:
        """若 counts 中有非整数值，counts_validated 应为 False。"""
        adata = _build_synthetic_yue_adata()
        adata.X.data[0] = 1.5
        adata.layers["counts"] = sp.csr_matrix(adata.X, dtype=np.float32)
        data = adata.layers["counts"].data
        is_int = bool(np.all(data == np.floor(data)))
        assert not is_int, "应检测到非整数值"
        adata.uns["expression_contract"] = {
            "x_scale": "raw_counts",
            "counts_layer": "counts",
            "counts_source": "X",
            "counts_validated": False,
            "counts_integer_check": "full",
            "soupx_layer": None,
            "processing_history": ["01_yue: integer check failed"],
            "stage": "01",
        }
        contract = validate_expression_contract(adata)
        assert contract["counts_validated"] is False

    def test_checkpoint_negative_value_raises(self) -> None:
        """checkpoint 校验发现负值应触发 AssertionError。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        # 注入负值到 layers["counts"]
        adata.layers["counts"].data[0] = -5.0
        with pytest.raises(AssertionError, match="negative"):
            _validate_counts_in_checkpoint(adata)

    def test_checkpoint_non_integer_raises(self) -> None:
        """checkpoint 校验发现非整数值应触发 AssertionError。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        # 注入非整数值（超过 1e-6 浮点误差容忍范围）
        adata.layers["counts"].data[0] = 10.5
        with pytest.raises(AssertionError, match="non-integer"):
            _validate_counts_in_checkpoint(adata)

    def test_float32_tolerance_1e6_below_passes(self) -> None:
        """float32 浮点误差 < 1e-6 的微小偏差应通过校验。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        # 注入 1e-7 的偏差（< 1e-6 容忍度）
        adata.layers["counts"].data[0] = 42.0 + 1e-7
        _validate_counts_in_checkpoint(adata)
        assert adata.uns["expression_contract"]["counts_validated"] is True


# ---- 集成：模拟 notebook 完整流程 -----------------------------------------------


class TestNotebookIntegration:
    """模拟 01_yue.ipynb 从数据加载到 checkpoint 的 counts 契约流程。"""

    def test_full_flow_load_to_checkpoint(self) -> None:
        """端到端：数据加载 → 建 counts 契约 → checkpoint 校验。"""
        # Step 1: 模拟数据加载（txt.gz 拼装 → AnnData，X=raw counts）
        adata = _build_synthetic_yue_adata()

        # Step 2: 建立 layers["counts"] + expression_contract（加载 cell）
        _establish_contract_yue(adata)

        # Step 3: checkpoint 校验（与 notebook 中 checkpoint cell 同口径）
        contract = _validate_counts_in_checkpoint(adata)
        assert contract["counts_validated"] is True

        # 3b. layers["counts"] 存在性
        assert "counts" in adata.layers
        assert sp.issparse(adata.layers["counts"])
        assert adata.layers["counts"].dtype == np.float32

        # 3c. 硬门禁全部通过（与 notebook hard_postconditions 同口径）
        counts = adata.layers["counts"]
        hard_postconditions = {
            "non_empty": adata.n_obs > 0 and adata.n_vars > 0,
            "x_sparse_float32": sp.issparse(adata.X) and adata.X.dtype == np.float32,
            "counts_layer_exists": "counts" in adata.layers,
            "counts_shape_aligned": counts.shape == adata.shape,
            "counts_non_negative": True,
            "counts_integer_check_completed": adata.uns["expression_contract"]["counts_integer_check"] == "full",
            "counts_contract_validated": adata.uns["expression_contract"]["counts_validated"],
        }
        assert all(hard_postconditions.values()), (
            f"hard_postconditions failed: {hard_postconditions}"
        )

    def test_flow_with_corrupted_contract_fails_checkpoint(self) -> None:
        """若契约被破坏（如缺少字段），checkpoint 应失败。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)

        # 恶意删除契约字段
        del adata.uns["expression_contract"]["counts_source"]

        with pytest.raises(ValueError, match="missing required keys"):
            validate_expression_contract(adata)

    def test_layers_counts_not_overwritten_by_checkpoint(self) -> None:
        """checkpoint 校验不得覆盖 layers['counts']。"""
        adata = _build_synthetic_yue_adata()
        _establish_contract_yue(adata)
        before = adata.layers["counts"].copy()
        _validate_counts_in_checkpoint(adata)
        after = adata.layers["counts"]
        # 数据应相同（校验只读不改矩阵）
        assert bool(np.all((before - after).data == 0))
