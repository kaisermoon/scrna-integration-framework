"""P1-d：02_merged.ipynb 的 doublet_include 过滤、每样本诊断与状态路由测试。

覆盖：
A. 算法参考实现锁定（_compute_doublet_inclusion）：
  1. 过滤正确性：include==False 被移除、True/NaN 保留
  2. 高 doublet 比例样本被 flag
  3. 高 uncertain 比例样本被 flag
  4. 小样本被 flag
  5. NaN include 样本被 flag
  6. doublet_prediction 缺失时降级
  7. 无异常路径：needs_review=False

B. 状态机接线（run_contract helpers）：
  8. needs_review=True -> StageStatus.NEEDS_REVIEW
  9. needs_review=False + all true -> StageStatus.SUCCESS
  10. hard_postcondition 某项 False -> FAILED（优先于 needs_review）
  11. promote_run 对 NEEDS_REVIEW manifest 抛错

C. Notebook 静态结构检查：
  12. 新 code cell 含 marker 与关键变量
  13. PARAMS cell 含 6 个 doublet 参数
  14. checkpoint cell 含 needs_review 路由 + doublet_inclusion
  15. 回归守卫：cell 0834b4b3 的 layers['counts'] 契约未被破坏
  16. 新 cell 不含 dir()
  17. notebook 不重定义 run_contract 函数
"""

from __future__ import annotations

import ast
import json
import tempfile
from pathlib import Path
from typing import Any

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scrna_integration.run_contract import (
    StageStatus,
    atomic_write_json,
    determine_stage_status,
    prepare_run,
    promote_run,
    sha256_file,
)

ROOT = Path(__file__).resolve().parents[1]
NB_PATH = ROOT / "notebooks" / "02_merged.ipynb"


# ---- 辅助函数 -------------------------------------------------------------------


def _load_nb_cells() -> list[dict[str, Any]]:
    with open(NB_PATH, encoding="utf-8") as f:
        return json.load(f)["cells"]


def _cell_source(cells: list[dict[str, Any]], cell_id: str) -> str:
    for c in cells:
        if c.get("id") == cell_id and c["cell_type"] == "code":
            return "".join(c["source"])
    raise ValueError(f"code cell id={cell_id} not found")


def _cell_by_marker(cells: list[dict[str, Any]], marker: str) -> str:
    """返回首个包含 marker 字符串的 code cell 源码。"""
    for c in cells:
        if c["cell_type"] == "code":
            src = "".join(c["source"])
            if marker in src:
                return src
    raise ValueError(f"No code cell contains marker: {marker}")


def _make_adata_with_doublet(
    n_cells: int = 100,
    n_genes: int = 50,
    *,
    seed: int = 42,
    doublet_include: np.ndarray | None = None,
    doublet_prediction: np.ndarray | None = None,
    sample_ids: np.ndarray | None = None,
    source_datasets: np.ndarray | None = None,
) -> anndata.AnnData:
    """构造带 doublet 列的合成 AnnData。"""
    rng = np.random.default_rng(seed)
    counts = sp.csr_matrix(rng.poisson(5, size=(n_cells, n_genes)).astype(np.float32))
    obs = {}
    if doublet_include is not None:
        obs["doublet_include"] = doublet_include
    if doublet_prediction is not None:
        obs["doublet_prediction"] = doublet_prediction
    if sample_ids is not None:
        obs["sample_id"] = sample_ids
    else:
        obs["sample_id"] = np.full(n_cells, "sample_A")
    if source_datasets is not None:
        obs["source_dataset"] = source_datasets
    else:
        obs["source_dataset"] = np.full(n_cells, "dataset_X")
    # 确保 obs DataFrame index 正确
    obs_df = pd.DataFrame(obs, index=[f"cell_{i}" for i in range(n_cells)])
    adata = anndata.AnnData(
        X=counts,
        obs=obs_df,
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)]),
        layers={"counts": counts.copy()},
        dtype=np.float32,
    )
    return adata


def _compute_doublet_inclusion(
    adata: anndata.AnnData,
    *,
    sample_key: str = "sample_id",
    include_key: str = "doublet_include",
    state_key: str = "doublet_prediction",
    max_doublet_frac: float = 0.30,
    max_uncertain_frac: float = 0.30,
    min_cells: int = 50,
) -> dict[str, Any]:
    """复刻 notebook 中 doublet 纳入构建的算法逻辑。

    此函数为算法真相参考，测试用例以此为准。
    返回 dict 含：filtered_adata, per_sample, flagged_samples, needs_review, report。
    """
    obs = adata.obs
    # 检查 include 列
    if include_key not in obs.columns:
        raise ValueError(f"obs missing {include_key}")

    # 确定分组键
    group_key = sample_key if sample_key in obs.columns else "source_dataset"

    # 每样本诊断
    per_sample = []
    for grp, grp_mask in obs.groupby(group_key, observed=False).groups.items():
        grp_obs = obs.loc[grp_mask]
        n_cells = len(grp_obs)
        include_mask = grp_obs[include_key]
        excluded_mask = include_mask == False  # noqa: E712
        n_excluded = int(excluded_mask.sum())
        n_included = n_cells - n_excluded
        doublet_frac = float(n_excluded / max(n_cells, 1))

        if state_key in grp_obs.columns:
            uncertain_mask = grp_obs[state_key] == "uncertain"
            uncertain_frac = float(uncertain_mask.sum() / max(n_cells, 1))
        else:
            uncertain_frac = np.nan

        flagged = False
        reasons = []
        if n_cells < min_cells:
            flagged = True
            reasons.append(f"too_few_cells({n_cells}<{min_cells})")
        if doublet_frac > max_doublet_frac:
            flagged = True
            reasons.append(f"doublet_frac({doublet_frac:.3f}>{max_doublet_frac})")
        if not np.isnan(uncertain_frac) and uncertain_frac > max_uncertain_frac:
            flagged = True
            reasons.append(f"uncertain_frac({uncertain_frac:.3f}>{max_uncertain_frac})")
        has_nan = include_mask.isna().any()
        if has_nan:
            flagged = True
            reasons.append("doublet_include_has_NaN(unstable_threshold)")

        src_ds = grp_obs["source_dataset"].iloc[0] if "source_dataset" in grp_obs.columns else str(grp)
        per_sample.append({
            "sample": str(grp),
            "source_dataset": str(src_ds),
            "n_cells": n_cells,
            "n_included": n_included,
            "n_excluded": n_excluded,
            "doublet_frac": round(doublet_frac, 4),
            "uncertain_frac": round(uncertain_frac, 4) if not np.isnan(uncertain_frac) else np.nan,
            "flagged": flagged,
            "flag_reason": "; ".join(reasons) if reasons else "",
        })

    # 过滤
    n_before = adata.n_obs
    keep_mask = (obs[include_key] != False)  # noqa: E712
    filtered = adata[keep_mask.values].copy()
    n_after = filtered.n_obs
    n_excluded_total = n_before - n_after

    # 保持 CSR + float32
    if sp.issparse(filtered.X):
        if filtered.X.dtype != np.float32:
            filtered.X = filtered.X.astype(np.float32)
    else:
        filtered.X = sp.csr_matrix(filtered.X, dtype=np.float32)

    flagged_samples = [d["sample"] for d in per_sample if d["flagged"]]
    needs_review = len(flagged_samples) > 0
    report = {
        "strategy": "exclude high-confidence doublets (doublet_include==False); keep singlet+uncertain",
        "include_key": include_key,
        "n_before": n_before,
        "n_after": n_after,
        "n_excluded": n_excluded_total,
        "per_sample": per_sample,
        "flagged_samples": flagged_samples,
        "needs_review": needs_review,
        "thresholds": {
            "max_doublet_fraction": max_doublet_frac,
            "max_uncertain_fraction": max_uncertain_frac,
            "min_sample_cells": min_cells,
        },
    }
    return {
        "filtered_adata": filtered,
        "per_sample": per_sample,
        "flagged_samples": flagged_samples,
        "needs_review": needs_review,
        "report": report,
    }


# ---- A. 算法参考实现测试 -----------------------------------------------------


class TestDoubletFiltering:
    """过滤正确性测试。"""

    def test_remove_false_keep_true_and_nan(self) -> None:
        """include==False 被移除；True 和 NaN 保留。"""
        include = np.array([True, True, False, False, True, np.nan, np.nan, False], dtype=object)
        adata = _make_adata_with_doublet(n_cells=8, doublet_include=include)
        result = _compute_doublet_inclusion(adata)
        filtered = result["filtered_adata"]
        # 5 cells kept: indices 0,1,4,5,6 (True + NaN)
        assert filtered.n_obs == 5, f"expected 5, got {filtered.n_obs}"
        assert result["report"]["n_before"] == 8
        assert result["report"]["n_after"] == 5
        assert result["report"]["n_excluded"] == 3
        # 验证保留的 cell 确实是 True 或 NaN
        kept_obs = filtered.obs["doublet_include"]
        assert kept_obs.isna().sum() == 2  # NaN 保留
        assert (kept_obs == True).sum() == 3  # noqa: E712

    def test_sparse_float32_preserved_after_filter(self) -> None:
        """过滤后 X 仍为 CSR + float32。"""
        include = np.array([True, False, True, True], dtype=object)
        adata = _make_adata_with_doublet(n_cells=4, doublet_include=include)
        result = _compute_doublet_inclusion(adata)
        filtered = result["filtered_adata"]
        assert sp.issparse(filtered.X), "X should remain sparse after filter"
        assert filtered.X.dtype == np.float32, f"X dtype {filtered.X.dtype} != float32"

    def test_counts_layer_preserved_after_filter(self) -> None:
        """过滤后 layers['counts'] 依然存在且 shape 与 X 一致。"""
        include = np.array([True, False, True, True, False], dtype=object)
        adata = _make_adata_with_doublet(n_cells=5, doublet_include=include)
        result = _compute_doublet_inclusion(adata)
        filtered = result["filtered_adata"]
        assert "counts" in filtered.layers, "counts layer should be preserved"
        assert filtered.layers["counts"].shape == filtered.X.shape


class TestFlagHighDoubletFraction:
    """高 doublet 比例样本被 flag。"""

    def test_flag_when_doublet_frac_exceeds_threshold(self) -> None:
        """doublet_frac > 0.30 -> flagged。"""
        # 40 cells, 15 False -> 0.375 > 0.30
        include = np.array([True] * 25 + [False] * 15, dtype=object)
        adata = _make_adata_with_doublet(n_cells=40, doublet_include=include)
        result = _compute_doublet_inclusion(adata, max_doublet_frac=0.30)
        per = result["per_sample"][0]
        assert per["flagged"] is True
        assert "doublet_frac" in per["flag_reason"]
        assert result["needs_review"] is True

    def test_no_flag_when_doublet_frac_within_threshold(self) -> None:
        """doublet_frac <= 0.30 -> not flagged。"""
        # 100 cells, 10 False -> 0.10 <= 0.30
        include = np.array([True] * 90 + [False] * 10, dtype=object)
        adata = _make_adata_with_doublet(n_cells=100, doublet_include=include)
        result = _compute_doublet_inclusion(adata, max_doublet_frac=0.30)
        per = result["per_sample"][0]
        # flagged=False unless other conditions trigger
        # n_cells=100 >= 50, doublet_frac=0.10 <= 0.30, no NaN
        assert per["flagged"] is False


class TestFlagHighUncertainFraction:
    """高 uncertain 比例样本被 flag。"""

    def test_flag_when_uncertain_frac_exceeds_threshold(self) -> None:
        """uncertain_frac > 0.30 -> flagged。"""
        include = np.full(60, True, dtype=object)
        prediction = np.array(["singlet"] * 20 + ["uncertain"] * 25 + ["doublet"] * 15)
        # 25/60 = 0.417 > 0.30
        adata = _make_adata_with_doublet(
            n_cells=60, doublet_include=include, doublet_prediction=prediction
        )
        result = _compute_doublet_inclusion(adata, max_uncertain_frac=0.30)
        per = result["per_sample"][0]
        assert per["flagged"] is True
        assert "uncertain_frac" in per["flag_reason"]


class TestFlagSmallSample:
    """小样本被 flag。"""

    def test_flag_when_fewer_than_min_cells(self) -> None:
        """n_cells < 50 -> flagged with too_few_cells reason。"""
        include = np.full(30, True, dtype=object)
        adata = _make_adata_with_doublet(n_cells=30, doublet_include=include)
        result = _compute_doublet_inclusion(adata, min_cells=50)
        per = result["per_sample"][0]
        assert per["flagged"] is True
        assert "too_few_cells" in per["flag_reason"]


class TestFlagNaNInclude:
    """NaN include 样本被 flag。"""

    def test_flag_when_has_nan_include(self) -> None:
        """样本含 doublet_include NaN -> flagged，但细胞仍保留。"""
        include = np.array([True] * 45 + [np.nan] * 5, dtype=object)
        adata = _make_adata_with_doublet(n_cells=50, doublet_include=include)
        result = _compute_doublet_inclusion(adata)
        per = result["per_sample"][0]
        assert per["flagged"] is True
        assert "doublet_include_has_NaN" in per["flag_reason"]
        # NaN 细胞仍保留
        filtered = result["filtered_adata"]
        assert filtered.n_obs == 50, "All cells should be kept (NaN preserved)"


class TestMissingPredictionColumn:
    """doublet_prediction 列缺失时降级。"""

    def test_missing_state_column_degraded(self) -> None:
        """不传 state 列 -> uncertain_frac=NaN、不报错、过滤仍生效。"""
        include = np.array([True] * 80 + [False] * 20, dtype=object)
        adata = _make_adata_with_doublet(n_cells=100, doublet_include=include)
        # 移除 doublet_prediction 列（构造时没传就是没有）
        # 我们构造一个没有 state 列的 adata
        # 需要先构造再删除
        rng = np.random.default_rng(42)
        counts = sp.csr_matrix(rng.poisson(5, size=(100, 50)).astype(np.float32))
        obs = pd.DataFrame({
            "doublet_include": include,
            "sample_id": "sample_A",
            "source_dataset": "dataset_X",
        }, index=[f"cell_{i}" for i in range(100)])
        adata_no_state = anndata.AnnData(
            X=counts, obs=obs,
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(50)]),
            dtype=np.float32,
        )
        result = _compute_doublet_inclusion(adata_no_state)
        per = result["per_sample"][0]
        assert np.isnan(per["uncertain_frac"]), "uncertain_frac should be NaN when column missing"
        # 过滤仍生效
        assert result["report"]["n_excluded"] == 20


class TestNoAbnormalSamples:
    """无异常路径。"""

    def test_all_normal_no_needs_review(self) -> None:
        """所有样本正常 -> needs_review=False。"""
        include = np.full(200, True, dtype=object)
        prediction = np.full(200, "singlet", dtype=object)
        adata = _make_adata_with_doublet(
            n_cells=200, doublet_include=include, doublet_prediction=prediction
        )
        result = _compute_doublet_inclusion(adata)
        assert result["needs_review"] is False
        assert result["flagged_samples"] == []


# ---- B. 状态机接线测试 -------------------------------------------------------


class TestStageStatusNeedsReview:
    """determine_stage_status 处理 needs_review 参数。"""

    def test_needs_review_true_yields_NEEDS_REVIEW(self) -> None:
        """needs_review=True 且 all true -> StageStatus.NEEDS_REVIEW。"""
        status = determine_stage_status(
            {}, {"a": True, "b": True},
            needs_review=True, allow_no_required_methods=True,
        )
        assert status == StageStatus.NEEDS_REVIEW

    def test_needs_review_false_all_true_yields_SUCCESS(self) -> None:
        """needs_review=False 且 all true -> StageStatus.SUCCESS。"""
        status = determine_stage_status(
            {}, {"a": True, "b": True},
            needs_review=False, allow_no_required_methods=True,
        )
        assert status == StageStatus.SUCCESS

    def test_failed_has_priority_over_needs_review(self) -> None:
        """hard_postcondition 某项 False -> FAILED，优先级高于 needs_review。"""
        status = determine_stage_status(
            {}, {"a": True, "b": False},
            needs_review=True, allow_no_required_methods=True,
        )
        assert status == StageStatus.FAILED


class TestPromoteRunBlocksNeedsReview:
    """promote_run 对 NEEDS_REVIEW manifest 抛错。"""

    def test_promote_run_raises_on_needs_review(self) -> None:
        """promote_run 对 NEEDS_REVIEW 状态抛 ValueError。"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "runs"
            run_id = "test-run-001"
            paths = prepare_run(root, run_id)

            # 写入 NEEDS_REVIEW manifest + dummy checkpoint
            checkpoint_path = paths.draft_dir / "dummy.h5ad"
            checkpoint_path.write_bytes(b"fake checkpoint")
            manifest = {
                "run_id": run_id,
                "stage": "02_merged",
                "stage_status": "NEEDS_REVIEW",
                "checkpoint": {
                    "path": "dummy.h5ad",
                    "sha256": sha256_file(checkpoint_path),
                },
            }
            atomic_write_json(paths.manifest_path, manifest)

            with pytest.raises(ValueError, match="cannot be promoted"):
                promote_run(paths)


# ---- C. Notebook 静态结构检查 -------------------------------------------------


class TestNotebookStructure:
    """Notebook 结构合规性检查。"""

    @pytest.fixture(scope="class")
    def cells(self) -> list[dict[str, Any]]:
        return _load_nb_cells()

    def test_new_code_cell_has_marker_and_vars(self, cells: list[dict[str, Any]]) -> None:
        """新 cell 含 marker '# === doublet 纳入构建' 与关键变量。"""
        src = _cell_by_marker(cells, "# === doublet 纳入构建 + 每样本异常诊断（决策8）===")
        assert "DOUBLET_INCLUDE_KEY" in src, "missing DOUBLET_INCLUDE_KEY reference"
        assert "SAMPLE_KEY" in src or "_sample_group_key" in src, "missing SAMPLE_KEY grouping"
        assert "_keep_mask" in src, "missing _keep_mask filter"
        assert "doublet_needs_review" in src, "missing doublet_needs_review"
        assert "flagged_samples" in src, "missing flagged_samples"

    def test_params_cell_has_doublet_params(self, cells: list[dict[str, Any]]) -> None:
        """PARAMS cell (id=9a522bbb) 含 6 个 doublet 参数。"""
        src = _cell_source(cells, "9a522bbb")
        required = [
            "DOUBLET_INCLUDE_KEY",
            "DOUBLET_STATE_KEY",
            "SAMPLE_KEY",
            "MAX_DOUBLET_FRACTION",
            "MAX_UNCERTAIN_FRACTION",
            "MIN_SAMPLE_CELLS_FOR_DOUBLET",
        ]
        for param in required:
            assert param in src, f"PARAMS cell missing {param}"

    def test_checkpoint_has_needs_review_routing(self, cells: list[dict[str, Any]]) -> None:
        """checkpoint cell (id=29c2b000) 含 needs_review 路由。"""
        src = _cell_source(cells, "29c2b000")
        assert "needs_review=doublet_needs_review" in src, (
            "determine_stage_status missing needs_review parameter"
        )
        assert "NEEDS_REVIEW" in src, "missing NEEDS_REVIEW string in checkpoint"
        assert "doublet_inclusion_report" in src, "missing doublet_inclusion_report"
        assert "doublet_include_present" in src, "hard_postconditions missing doublet_include_present"
        assert "doublet_filter_applied" in src, "hard_postconditions missing doublet_filter_applied"
        assert '"doublet_inclusion"' in src, (
            "manifest_payload/uns missing doublet_inclusion key"
        )

    def test_checkpoint_has_promote_run_only_in_else(self, cells: list[dict[str, Any]]) -> None:
        """promote_run 仅在 else 分支中调用，NEEDS_REVIEW 分支不调用。"""
        src = _cell_source(cells, "29c2b000")
        # 找到 NEEDS_REVIEW 分支后的代码并确认 promote_run 只在 else 里
        needs_idx = src.find('if stage_status.value == "NEEDS_REVIEW":')
        assert needs_idx > 0, "NEEDS_REVIEW branch not found"
        # NEEDS_REVIEW 分支内不应有 promote_run 调用（注释里的不算）
        # 找到下一个 else 的位置
        else_idx = src.find("else:", needs_idx)
        assert else_idx > needs_idx, "else branch not found after NEEDS_REVIEW"
        needs_block = src[needs_idx:else_idx]
        # 注释里可以有 promote_run 字样（"raise"说明），但不应有 promote_run( 的调用
        assert "promote_run(" not in needs_block, (
            "promote_run() should NOT be called in NEEDS_REVIEW branch"
        )
        # else 分支应该调用 promote_run
        else_block = src[else_idx:]
        assert "promote_run(" in else_block, (
            "promote_run() should be called in else (SUCCESS) branch"
        )

    def test_redline_cell_0834b4b3_untouched(self, cells: list[dict[str, Any]]) -> None:
        """回归守卫：cell 0834b4b3 仍含 layers['counts'] 契约 assert。"""
        src = _cell_source(cells, "0834b4b3")
        assert "layers[\"counts\"]" in src or "layers['counts']" in src, (
            "REDLINE: cell 0834b4b3 no longer references layers['counts']"
        )
        assert "non-negative integers" in src, (
            "REDLINE: cell 0834b4b3 counts_integer_check removed"
        )

    def test_redline_expression_contract_preserved(self, cells: list[dict[str, Any]]) -> None:
        """回归守卫：checkpoint cell 仍含 expression_contract 写入 + validate。"""
        src = _cell_source(cells, "29c2b000")
        assert "# === expression_contract" in src, (
            "REDLINE: expression_contract comment removed from checkpoint"
        )
        assert "validate_expression_contract(adata" in src, (
            "REDLINE: validate_expression_contract call removed from checkpoint"
        )

    def test_no_dir_in_any_code_cell(self, cells: list[dict[str, Any]]) -> None:
        """所有 code cell 不含 dir() 调用。"""
        for c in cells:
            if c["cell_type"] == "code":
                src = "".join(c["source"])
                assert "dir()" not in src, (
                    f"dir() found in code cell (id={c.get('id', 'N/A')})"
                )

    def test_no_redefinition_of_run_contract_functions(self, cells: list[dict[str, Any]]) -> None:
        """notebook 不重定义 run_contract 已有函数。"""
        protected = [
            "determine_stage_status",
            "validate_expression_contract",
            "promote_run",
            "prepare_run",
        ]
        for func_name in protected:
            for c in cells:
                if c["cell_type"] == "code":
                    src = "".join(c["source"])
                    assert f"def {func_name}" not in src, (
                        f"notebook redefines {func_name} (should import from run_contract)"
                    )

    def test_all_code_cells_parse(self, cells: list[dict[str, Any]]) -> None:
        """每个 code cell 独立通过 ast.parse。"""
        for c in cells:
            if c["cell_type"] == "code":
                src = "".join(c["source"])
                try:
                    ast.parse(src)
                except SyntaxError as e:
                    pytest.fail(f"SyntaxError in cell id={c.get('id', 'N/A')}: {e}")

    def test_notebook_json_valid(self) -> None:
        """notebook JSON 合法，有 cells 数组。"""
        assert NB_PATH.exists(), f"not found: {NB_PATH}"
        with open(NB_PATH, encoding="utf-8") as f:
            nb = json.load(f)
        assert "cells" in nb, "notebook missing 'cells' key"
        assert len(nb["cells"]) >= 31, f"expected >=31 cells, got {len(nb['cells'])}"

    def test_new_markdown_cell_present(self, cells: list[dict[str, Any]]) -> None:
        """新 markdown cell 含有关于 doublet_include 构建的说明。"""
        md_found = False
        for c in cells:
            if c["cell_type"] == "markdown":
                src = "".join(c["source"])
                if "doublet_include" in src and "构建正式整合对象" in src:
                    md_found = True
                    break
        assert md_found, "new markdown cell about doublet inclusion not found"
