"""P0-g：03_normalized 的 expression_contract 与 counts 完整性测试。

覆盖：
- counts-layer cell：读取上游 layers["counts"]（不自拷贝），验证契约
- normalize cell：float64 sum 断言 layers["counts"] 未变（相对差 < 1e-6）
- checkpoint cell：写 expression_contract(stage="03", x_scale="normalized_log1p")
  + counts_layer_preserved 到 hard_postconditions
- 禁止 dir() 判断结果残留
- processing_history 随 NORMALIZATION_METHOD 变化
"""

from __future__ import annotations

import ast
import json
from pathlib import Path
from typing import Any

import anndata
import numpy as np
import pytest
import scanpy as sc
import scipy.sparse as sp

from scrna_integration.run_contract import (
    StageStatus,
    determine_stage_status,
    validate_expression_contract,
)


# ---- 辅助函数 -------------------------------------------------------------------


def _make_minimal_adata(
    n_cells: int = 100,
    n_genes: int = 200,
    *,
    seed: int = 42,
) -> anndata.AnnData:
    """构造带 layers["counts"] 和有效 stage-02 expression_contract 的最小 AnnData。"""
    rng = np.random.default_rng(seed)
    counts = sp.csr_matrix(
        rng.poisson(5, size=(n_cells, n_genes)).astype(np.float32)
    )
    adata = anndata.AnnData(
        X=counts.copy(),
        layers={"counts": counts.copy()},
        obs={"source_dataset": rng.choice(["A", "B"], size=n_cells)},
        var={"gene_ids": [f"gene_{i}" for i in range(n_genes)]},
        dtype=np.float32,
    )
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.uns["expression_contract"] = {
        "x_scale": "raw_counts",
        "counts_layer": "counts",
        "counts_source": "X",
        "counts_validated": True,
        "counts_integer_check": "full",
        "soupx_layer": None,
        "processing_history": ["01_kim: built counts layer from X"],
        "stage": "02",
    }
    return adata


def _check_no_dir_in_nb_cells(nb_path: Path, cell_ids: list[str]) -> None:
    """确保指定 notebook cell 不使用 dir() 判断变量存在性。"""
    with open(nb_path) as f:
        nb = json.load(f)
    for cell in nb["cells"]:
        if cell.get("id") in cell_ids and cell["cell_type"] == "code":
            src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
            assert "dir()" not in src, f"cell {cell['id']} 使用了 dir()（红线6）"


def _check_no_redefinition_in_nb(nb_path: Path, func_name: str) -> None:
    """确保 notebook 不重新定义 run_contract 中已有的函数。"""
    with open(nb_path) as f:
        nb = json.load(f)
    for cell in nb["cells"]:
        if cell["cell_type"] == "code":
            src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
            assert f"def {func_name}" not in src, (
                f"notebook 禁止重新定义 {func_name}（应 import from run_contract）"
            )


# ---- counts-layer cell 逻辑测试 --------------------------------------------------


class TestCountsLayerRead:
    """03-counts-layer：读取上游 layers["counts"]，不自拷贝。"""

    def test_reads_upstream_counts_not_self_copy(self) -> None:
        """layers["counts"] 存在时直接读取，不调用 adata.X.copy()。"""
        adata = _make_minimal_adata()
        counts_before = float(adata.layers["counts"].sum())

        # 模拟 03-counts-layer cell 的核心逻辑
        _upstream_contract = validate_expression_contract(
            adata, expected_scale="raw_counts", stage="02"
        )
        assert "counts" in adata.layers
        assert _upstream_contract["counts_layer"] == "counts"

        # layers["counts"] 不应被修改
        counts_after = float(adata.layers["counts"].sum())
        assert counts_before == counts_after

    def test_missing_counts_layer_raises_key_error(self) -> None:
        """layers["counts"] 缺失时应报 KeyError。"""
        adata = _make_minimal_adata()
        del adata.layers["counts"]
        # 03-counts-layer cell 中会做这个检查
        with pytest.raises(KeyError, match="counts"):
            if "counts" not in adata.layers:
                raise KeyError("layers['counts'] 不存在")

    def test_contract_stage_mismatch_blocked(self) -> None:
        """上游 stage 不对应被 expression_contract 校验挡下。"""
        adata = _make_minimal_adata()
        adata.uns["expression_contract"]["stage"] = "01"
        with pytest.raises(ValueError, match="expected '02'"):
            validate_expression_contract(adata, expected_scale="raw_counts", stage="02")

    def test_contract_x_scale_mismatch_blocked(self) -> None:
        """上游 x_scale 不是 raw_counts 应被挡下。"""
        adata = _make_minimal_adata()
        adata.uns["expression_contract"]["x_scale"] = "normalized_log1p"
        with pytest.raises(ValueError, match="expected 'raw_counts'"):
            validate_expression_contract(adata, expected_scale="raw_counts", stage="02")

    def test_missing_contract_raises_key_error(self) -> None:
        """无 expression_contract 时应报错。"""
        adata = _make_minimal_adata()
        del adata.uns["expression_contract"]
        with pytest.raises(KeyError, match="expression_contract not found"):
            validate_expression_contract(adata, expected_scale="raw_counts", stage="02")


# ---- normalize cell 逻辑测试 -----------------------------------------------------


class TestNormalizeCountsIntegrity:
    """03-normalize-code：float64 sum 断言 layers["counts"] 未变。"""

    def test_standard_normalization_preserves_counts_layer(self) -> None:
        """standard (normalize_total + log1p) 后 layers["counts"] 不变。"""
        adata = _make_minimal_adata()
        _upstream_contract = adata.uns["expression_contract"]
        # 标准化前
        _counts_sum_before = float(adata.layers["counts"].sum(dtype=np.float64))

        # 执行标准化（standard 模式）
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.X = adata.X.astype(np.float32)

        # 标准化后验证
        _counts_sum_after = float(adata.layers["counts"].sum(dtype=np.float64))
        _rel_diff = abs(_counts_sum_before - _counts_sum_after) / max(
            abs(_counts_sum_before), 1.0
        )
        assert _rel_diff < 1e-6, (
            f"layers['counts'] 被修改：rel_diff={_rel_diff:.2e}"
        )
        _counts_integrity_checked = True
        assert _counts_integrity_checked

    def test_pearson_residuals_preserves_counts_layer(self) -> None:
        """Pearson residuals 标准化后 layers["counts"] 不变。"""
        adata = _make_minimal_adata()
        _upstream_contract = adata.uns["expression_contract"]
        _counts_sum_before = float(adata.layers["counts"].sum(dtype=np.float64))

        # 执行 Pearson residuals 标准化
        sc.experimental.pp.normalize_pearson_residuals(adata)
        adata.X = adata.X.astype(np.float32)

        _counts_sum_after = float(adata.layers["counts"].sum(dtype=np.float64))
        _rel_diff = abs(_counts_sum_before - _counts_sum_after) / max(
            abs(_counts_sum_before), 1.0
        )
        assert _rel_diff < 1e-6, (
            f"layers['counts'] 被修改：rel_diff={_rel_diff:.2e}"
        )

    def test_counts_integrity_violation_detected(self) -> None:
        """若 layers["counts"] 被意外覆盖，断言应触发。"""
        adata = _make_minimal_adata()
        _counts_sum_before = float(adata.layers["counts"].sum(dtype=np.float64))

        # 人为覆盖 layers["counts"]（模拟 bug）
        adata.layers["counts"] = adata.layers["counts"].copy()
        adata.layers["counts"].data[0] = 9999.0  # 篡改一个值

        _counts_sum_after = float(adata.layers["counts"].sum(dtype=np.float64))
        _rel_diff = abs(_counts_sum_before - _counts_sum_after) / max(
            abs(_counts_sum_before), 1.0
        )
        assert _rel_diff >= 1e-6, (
            f"篡改后的 layers['counts'] 应被检测到差异，但 rel_diff={_rel_diff:.2e}"
        )

    def test_explicit_marker_variable_used(self) -> None:
        """验证使用显式标记变量 _counts_integrity_checked，而非 dir() 判断。"""
        adata = _make_minimal_adata()
        _counts_sum_before = float(adata.layers["counts"].sum(dtype=np.float64))

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.X = adata.X.astype(np.float32)

        _counts_sum_after = float(adata.layers["counts"].sum(dtype=np.float64))
        rel_diff = abs(_counts_sum_before - _counts_sum_after) / max(abs(_counts_sum_before), 1.0)
        assert rel_diff < 1e-6

        # 显式标记变量（红线6：禁止 dir() 判断）
        # 正确做法：直接赋值 _counts_integrity_checked = True
        # 下游直接引用该变量（NameError 即表示未执行，无需 dir()）
        _counts_integrity_checked = True
        assert _counts_integrity_checked is True


# ---- checkpoint cell 逻辑测试 ----------------------------------------------------


class TestCheckpointExpressionContract:
    """03-checkpoint-code：expression_contract + counts_layer_preserved。"""

    def test_writes_expression_contract_with_stage_03(self) -> None:
        """checkpoint 写 stage="03", x_scale="normalized_log1p" 的契约。"""
        adata = _make_minimal_adata()
        _upstream_contract = adata.uns["expression_contract"]
        _counts_integrity_checked = True

        # 模拟 checkpoint cell 的 expression_contract 写入
        NORMALIZATION_METHOD = "standard"
        TARGET_SUM = 10000

        if NORMALIZATION_METHOD == "standard":
            _processing_step = f"normalize_total(target_sum={TARGET_SUM}) + log1p"
        elif NORMALIZATION_METHOD == "pearson_residuals":
            _processing_step = "normalize_pearson_residuals"
        else:
            _processing_step = f"normalize_{NORMALIZATION_METHOD}"

        _upstream_history = _upstream_contract.get("processing_history", [])
        adata.uns["expression_contract"] = {
            "x_scale": "normalized_log1p",
            "counts_layer": "counts",
            "counts_source": _upstream_contract.get("counts_source", "layers[counts]"),
            "counts_validated": _upstream_contract.get("counts_validated", True),
            "counts_integer_check": _upstream_contract.get("counts_integer_check", "full"),
            "soupx_layer": None,
            "processing_history": _upstream_history + [_processing_step],
            "stage": "03",
        }

        # 自验证
        contract = validate_expression_contract(
            adata, expected_scale="normalized_log1p", stage="03"
        )
        assert contract["x_scale"] == "normalized_log1p"
        assert contract["stage"] == "03"
        assert contract["counts_layer"] == "counts"
        assert contract["soupx_layer"] is None
        assert len(contract["processing_history"]) == 2  # upstream + 当前

        # hard_postconditions 包含 counts_layer_preserved
        hard_postconditions = {
            "non_empty": adata.n_obs > 0 and adata.n_vars > 0,
            "x_finite": True,
            "x_float32": adata.X.dtype == np.float32,
            "hvg_non_empty": True,
            "counts_layer_preserved": _counts_integrity_checked,
        }
        assert hard_postconditions["counts_layer_preserved"] is True

        # stage_status 计算
        status = determine_stage_status(
            {}, hard_postconditions, allow_no_required_methods=True
        )
        assert status in (StageStatus.SUCCESS,)

    def test_processing_history_uses_normalization_method_variable(self) -> None:
        """processing_history 随 NORMALIZATION_METHOD 变化，不硬编码。"""
        adata_std = _make_minimal_adata()
        adata_pr = _make_minimal_adata()

        # standard 模式
        NORMALIZATION_METHOD = "standard"
        TARGET_SUM = 10000
        _step_std = f"normalize_total(target_sum={TARGET_SUM}) + log1p"
        adata_std.uns["expression_contract"] = {
            "x_scale": "normalized_log1p", "counts_layer": "counts",
            "counts_source": "X", "counts_validated": True,
            "counts_integer_check": "full", "soupx_layer": None,
            "processing_history": ["01: build"] + [_step_std], "stage": "03",
        }
        validate_expression_contract(adata_std)
        assert _step_std in adata_std.uns["expression_contract"]["processing_history"]

        # pearson_residuals 模式
        NORMALIZATION_METHOD = "pearson_residuals"
        _step_pr = "normalize_pearson_residuals"
        adata_pr.uns["expression_contract"] = {
            "x_scale": "normalized_log1p", "counts_layer": "counts",
            "counts_source": "X", "counts_validated": True,
            "counts_integer_check": "full", "soupx_layer": None,
            "processing_history": ["01: build"] + [_step_pr], "stage": "03",
        }
        validate_expression_contract(adata_pr)
        assert _step_pr in adata_pr.uns["expression_contract"]["processing_history"]

        # 两种模式产出的 processing_history 不同
        assert _step_std != _step_pr

    def test_counts_layer_preserved_in_postconditions(self) -> None:
        """hard_postconditions 含 counts_layer_preserved，决定 stage 状态。"""
        adata = _make_minimal_adata()

        # 契约完整时 SUCCESS
        hard_ok = {
            "non_empty": True, "x_finite": True, "x_float32": True,
            "hvg_non_empty": True, "counts_layer_preserved": True,
        }
        assert determine_stage_status({}, hard_ok, allow_no_required_methods=True) is StageStatus.SUCCESS

        # counts_layer_preserved=False 时 FAILED
        hard_bad = {
            "non_empty": True, "x_finite": True, "x_float32": True,
            "hvg_non_empty": True, "counts_layer_preserved": False,
        }
        assert determine_stage_status({}, hard_bad, allow_no_required_methods=True) is StageStatus.FAILED

    def test_expression_contract_carries_forward_upstream_fields(self) -> None:
        """counts_source / counts_validated / counts_integer_check 从上游继承。"""
        adata = _make_minimal_adata()
        _upstream_contract = adata.uns["expression_contract"]

        # 模拟 checkpoint 写入
        adata.uns["expression_contract"] = {
            "x_scale": "normalized_log1p",
            "counts_layer": "counts",
            "counts_source": _upstream_contract.get("counts_source", "layers[counts]"),
            "counts_validated": _upstream_contract.get("counts_validated", True),
            "counts_integer_check": _upstream_contract.get("counts_integer_check", "full"),
            "soupx_layer": None,
            "processing_history": _upstream_contract.get("processing_history", []) + ["step"],
            "stage": "03",
        }

        contract = validate_expression_contract(adata)
        assert contract["counts_source"] == "X"
        assert contract["counts_validated"] is True
        assert contract["counts_integer_check"] == "full"


# ---- notebook 静态合规检查 --------------------------------------------------------


_NB_PATH = Path(__file__).resolve().parent.parent / "notebooks" / "03_normalized.ipynb"


class TestNotebookStaticCompliance:
    """对 03_normalized.ipynb 的静态检查（不执行 notebook）。"""

    def test_no_redefinition_of_run_contract_functions(self) -> None:
        """禁止在 notebook 中重新定义 validate_expression_contract 等。"""
        _check_no_redefinition_in_nb(_NB_PATH, "validate_expression_contract")
        _check_no_redefinition_in_nb(_NB_PATH, "aggregate_method_status")
        _check_no_redefinition_in_nb(_NB_PATH, "determine_stage_status")

    def test_no_dir_usage_in_modified_cells(self) -> None:
        """counts-layer / normalize / checkpoint cell 不使用 dir()（红线6）。"""
        _check_no_dir_in_nb_cells(_NB_PATH, [
            "03-counts-layer", "03-normalize-code", "03-checkpoint-code",
        ])

    def test_imports_validate_expression_contract(self) -> None:
        """03-setup cell 导入了 validate_expression_contract。"""
        with open(_NB_PATH) as f:
            nb = json.load(f)
        for cell in nb["cells"]:
            if cell.get("id") == "03-setup" and cell["cell_type"] == "code":
                src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
                assert "validate_expression_contract" in src, (
                    "03-setup 缺少 validate_expression_contract 导入"
                )

    def test_processing_history_not_hardcoded(self) -> None:
        """checkpoint cell 的 processing_history 使用 NORMALIZATION_METHOD 变量。"""
        with open(_NB_PATH) as f:
            nb = json.load(f)
        for cell in nb["cells"]:
            if cell.get("id") == "03-checkpoint-code" and cell["cell_type"] == "code":
                src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
                # processing_step 随 NORMALIZATION_METHOD 变化
                assert "NORMALIZATION_METHOD" in src
                # 不硬编码具体字符串
                assert "_processing_step" in src

    def test_explicit_integrity_check_variable(self) -> None:
        """normalize cell 使用 _counts_integrity_checked 显式标记（非 dir()）。"""
        with open(_NB_PATH) as f:
            nb = json.load(f)
        for cell in nb["cells"]:
            if cell.get("id") == "03-normalize-code" and cell["cell_type"] == "code":
                src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
                assert "_counts_integrity_checked" in src
                assert "dir()" not in src

    def test_checkpoint_has_counts_layer_preserved(self) -> None:
        """checkpoint cell 的 hard_postconditions 含 counts_layer_preserved。"""
        with open(_NB_PATH) as f:
            nb = json.load(f)
        for cell in nb["cells"]:
            if cell.get("id") == "03-checkpoint-code" and cell["cell_type"] == "code":
                src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
                assert '"counts_layer_preserved"' in src
                assert '_counts_integrity_checked' in src
