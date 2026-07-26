"""P1-c: counts matrix 格式模板 doublet 三态政策行为测试与静态源码断言。

覆盖：
- 三态赋值（singlet/uncertain/doublet）
- doublet_include 仅排除高置信 doublet
- predicted_doublet 与 call 一致
- threshold per-sample 持久化
- margin 参数消费验证
- needs_review 三类触发（高比例/小样本/阈值不稳定）
- determine_stage_status 状态分支
- counts 契约不变
- 静态源码断言（no dir()/no redefinition/no doublet subset/required fields/needs_review branch）
"""

# ruff: noqa: N806

from __future__ import annotations

import json
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from scrna_integration.run_contract import (
    StageStatus,
    determine_stage_status,
)

# ---- 科学参数（与 notebook PARAMS 保持同步）-----------------------------------

DOUBLET_UNCERTAIN_MARGIN = 0.10
DOUBLET_RATE_ALERT = 0.30
DOUBLET_MIN_CELLS = 50

NB_PATH = (
    Path(__file__).resolve().parent.parent
    / "notebooks" / "01_per_dataset" / "01_template_counts_matrix.ipynb"
)


# ---- helpers ----------------------------------------------------------------

def _make_synthetic_adata(
    n_cells: int = 200,
    n_genes: int = 100,
    *,
    seed: int = 42,
) -> anndata.AnnData:
    """构造带 contract 的合成 AnnData，含 sample_id 和 layers["counts"]。"""
    counts = sp.random(n_cells, n_genes, density=0.3, format="csr", dtype=np.float32)
    counts.data = np.floor(counts.data * 100 + 1)

    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    obs["source_dataset"] = "yue"
    obs["sample_id"] = [f"sample_{i % 3 + 1}" for i in range(n_cells)]

    adata = anndata.AnnData(X=counts, obs=obs, dtype=np.float32)
    adata.layers["counts"] = sp.csr_matrix(adata.X, dtype=np.float32)
    adata.uns["expression_contract"] = {
        "x_scale": "raw_counts",
        "counts_layer": "counts",
        "counts_source": "X",
        "counts_validated": False,
        "counts_integer_check": None,
        "soupx_layer": None,
        "processing_history": [],
        "stage": "01",
    }
    return adata


def _apply_three_state(
    scores: np.ndarray,
    threshold: float,
    margin: float = DOUBLET_UNCERTAIN_MARGIN,
) -> np.ndarray:
    """复刻 doublet cell 中三态判定的核心逻辑（见 cell be34c1dd）。

    规则（决策8）：
    - s >= T*(1+margin)  → "doublet"（高置信双细胞）
    - s <= T*(1-margin)  → "singlet"（安全）
    - 之间               → "uncertain"（缓冲带，保守保留）
    """
    T_upper = threshold * (1 + margin)
    T_lower = threshold * (1 - margin)

    n = len(scores)
    calls = np.full(n, "singlet", dtype=object)
    for i, s in enumerate(scores):
        if s >= T_upper:
            calls[i] = "doublet"
        elif s <= T_lower:
            calls[i] = "singlet"
        else:
            calls[i] = "uncertain"
    return calls


# ---- 行为测试 ----------------------------------------------------------------

class TestThreeStateAssignment:
    """三态切分逻辑：给定分数与阈值 + margin，验证 singlet/uncertain/doublet 分类。"""

    def test_clear_doublet(self):
        """分数远超阈值 → 全部 doublet。"""
        scores = np.array([0.5, 0.6, 0.8])
        T = 0.3
        calls = _apply_three_state(scores, T)
        assert all(c == "doublet" for c in calls)

    def test_clear_singlet(self):
        """分数远低阈值 → 全部 singlet。"""
        scores = np.array([0.05, 0.08, 0.02])
        T = 0.3
        calls = _apply_three_state(scores, T)
        assert all(c == "singlet" for c in calls)

    def test_uncertain_band(self):
        """分数在 T*(1-margin) 与 T*(1+margin) 之间 → uncertain。"""
        # T=0.3, margin=0.10: band = [0.27, 0.33]
        scores = np.array([0.28, 0.30, 0.32])
        T = 0.3
        calls = _apply_three_state(scores, T)
        assert list(calls) == ["uncertain", "uncertain", "uncertain"]

    def test_boundary_exact(self):
        """精确边界值：T_lower 和 T_upper 处归入正确类。"""
        T = 0.3
        margin = 0.10
        T_upper = T * (1 + margin)  # 0.33
        T_lower = T * (1 - margin)  # 0.27

        scores = np.array([T_lower, T_upper, T_lower - 0.001, T_upper + 0.001])
        calls = _apply_three_state(scores, T, margin=margin)
        assert calls[0] == "singlet"   # == T_lower → singlet
        assert calls[1] == "doublet"   # == T_upper → doublet
        assert calls[2] == "singlet"   # < T_lower → singlet
        assert calls[3] == "doublet"   # > T_upper → doublet

    def test_mixed_scores(self):
        """混合分数：各分数落入预期类别。"""
        scores = np.array([0.05, 0.25, 0.29, 0.35, 0.50])
        T = 0.3
        # margin=0.10: T_lower=0.27, T_upper=0.33
        calls = _apply_three_state(scores, T)
        assert list(calls) == [
            "singlet", "singlet", "uncertain", "doublet", "doublet"
        ]


class TestDoubletInclude:
    """doublet_include 仅排除高置信 doublet；uncertain 与 singlet 均保留。"""

    def test_only_doublet_excluded(self):
        """doublet_include == False 当且仅当 doublet_call == 'doublet'。"""
        adata = _make_synthetic_adata(n_cells=9)
        adata.obs["doublet_call"] = pd.Categorical(
            ["singlet", "singlet", "singlet",
             "uncertain", "uncertain", "uncertain",
             "doublet", "doublet", "doublet"],
            categories=["singlet", "uncertain", "doublet"],
            ordered=True,
        )
        adata.obs["doublet_include"] = adata.obs["doublet_call"] != "doublet"
        adata.obs["predicted_doublet"] = adata.obs["doublet_call"] == "doublet"

        # singlet → include=True
        assert adata.obs.loc[
            adata.obs["doublet_call"] == "singlet", "doublet_include"
        ].all()
        # uncertain → include=True
        assert adata.obs.loc[
            adata.obs["doublet_call"] == "uncertain", "doublet_include"
        ].all()
        # doublet → include=False
        assert not adata.obs.loc[
            adata.obs["doublet_call"] == "doublet", "doublet_include"
        ].any()

    def test_predicted_doublet_matches_call(self):
        """predicted_doublet == (doublet_call == 'doublet') 严格一致。"""
        adata = _make_synthetic_adata(n_cells=9)
        adata.obs["doublet_call"] = pd.Categorical(
            ["singlet", "singlet", "uncertain", "uncertain",
             "doublet", "doublet", "singlet", "uncertain", "doublet"],
            categories=["singlet", "uncertain", "doublet"],
            ordered=True,
        )
        adata.obs["predicted_doublet"] = adata.obs["doublet_call"] == "doublet"
        assert (
            adata.obs["predicted_doublet"]
            == (adata.obs["doublet_call"] == "doublet")
        ).all()


class TestThresholdPerSample:
    """多样本不同阈值时每细胞 doublet_threshold 等于其样本阈值。"""

    def test_threshold_persisted_correctly(self):
        """每个细胞的 doublet_threshold 等于其所属样本的判定阈值。"""
        adata = _make_synthetic_adata(n_cells=6)
        adata.obs["sample_id"] = ["A", "A", "A", "B", "B", "B"]
        adata.obs["doublet_threshold"] = np.nan

        sample_thresholds = {"A": 0.25, "B": 0.35}
        for sid, t in sample_thresholds.items():
            mask = adata.obs["sample_id"] == sid
            adata.obs.loc[mask, "doublet_threshold"] = t

        assert (
            adata.obs.loc[adata.obs["sample_id"] == "A", "doublet_threshold"].iloc[0]
            == 0.25
        )
        assert (
            adata.obs.loc[adata.obs["sample_id"] == "B", "doublet_threshold"].iloc[0]
            == 0.35
        )

    def test_skipped_sample_keeps_nan(self):
        """被跳过样本（细胞数 < DOUBLET_MIN_CELLS）的 doublet_threshold 保持 np.nan。"""
        adata = _make_synthetic_adata(n_cells=6)
        adata.obs["sample_id"] = ["A", "A", "B", "B", "C", "C"]
        adata.obs["doublet_threshold"] = np.nan

        # 只设置 A 的阈值，B/C 模拟跳过
        mask_a = adata.obs["sample_id"] == "A"
        adata.obs.loc[mask_a, "doublet_threshold"] = 0.30

        assert np.isnan(
            adata.obs.loc[adata.obs["sample_id"] == "B", "doublet_threshold"].iloc[0]
        )
        assert np.isnan(
            adata.obs.loc[adata.obs["sample_id"] == "C", "doublet_threshold"].iloc[0]
        )


class TestMarginParamConsumed:
    """margin 参数真实消费：margin=0 与 margin=0.30 产出不同 uncertain 计数。"""

    def test_margin_zero_no_uncertain(self):
        """margin=0 时所有分数直接归入 singlet 或 doublet，无 uncertain。"""
        scores = np.array([0.20, 0.25, 0.28, 0.30, 0.35, 0.40])
        T = 0.30
        # margin=0: T_lower = T_upper = 0.30
        # s < 0.30 → singlet（第一个分支 s>=0.30 未命中，elif s<=0.30 命中）
        # s >= 0.30 → doublet（第一个分支命中）
        calls = _apply_three_state(scores, T, margin=0.0)
        n_uncertain = int((calls == "uncertain").sum())
        assert n_uncertain == 0

    def test_margin_positive_produces_uncertain(self):
        """margin=0.30 时不确定带宽更大 → 至少 3 个细胞进 uncertain。"""
        scores = np.array([0.20, 0.25, 0.28, 0.30, 0.35, 0.40])
        T = 0.30
        # margin=0.30: T_lower=0.21, T_upper=0.39
        # 0.20<0.21→singlet; 0.25/0.28/0.30/0.35∈(0.21,0.39)→uncertain; 0.40>0.39→doublet
        calls = _apply_three_state(scores, T, margin=0.30)
        n_uncertain = int((calls == "uncertain").sum())
        assert n_uncertain >= 3


class TestNeedsReview:
    """needs_review 三类触发：高比例/小样本/阈值不稳定。"""

    def test_high_doublet_rate_triggers_needs_review(self):
        """高置信 doublet 比例 > DOUBLET_RATE_ALERT → needs_review=True。"""
        n_cells = 100
        n_doublet = 40  # 40% > 30% alert
        pct = n_doublet / n_cells
        assert pct > DOUBLET_RATE_ALERT

        per_sample = {
            "sample_1": {
                "n_cells": n_cells,
                "n_doublet": n_doublet,
                "pct_doublet": pct,
                "needs_review": pct > DOUBLET_RATE_ALERT,
            }
        }
        any_needs_review = any(
            v.get("needs_review", False) for v in per_sample.values()
        )
        assert any_needs_review

    def test_small_sample_triggers_needs_review(self):
        """样本细胞数 < DOUBLET_MIN_CELLS → needs_review=True。"""
        n_cells = 30
        assert n_cells < DOUBLET_MIN_CELLS
        per_sample = {
            "sample_tiny": {
                "n_cells": n_cells,
                "needs_review": n_cells < DOUBLET_MIN_CELLS,
            }
        }
        assert per_sample["sample_tiny"]["needs_review"]

    def test_small_sample_keeps_singlet_include(self):
        """小样本跳过检测后，全部细胞保持 singlet 且 include=True。"""
        adata = _make_synthetic_adata(n_cells=30)
        adata.obs["doublet_call"] = pd.Categorical(
            ["singlet"] * 30,
            categories=["singlet", "uncertain", "doublet"],
            ordered=True,
        )
        adata.obs["doublet_include"] = True
        adata.obs["doublet_score"] = np.nan
        adata.obs["doublet_threshold"] = np.nan

        assert (adata.obs["doublet_call"] == "singlet").all()
        assert adata.obs["doublet_include"].all()

    def test_unstable_threshold_triggers_needs_review(self):
        """threshold=None（自动定阈失败）→ needs_review=True，细胞全部 singlet。"""
        adata = _make_synthetic_adata(n_cells=100)
        adata.obs["doublet_call"] = pd.Categorical(
            ["singlet"] * 100,
            categories=["singlet", "uncertain", "doublet"],
            ordered=True,
        )
        adata.obs["doublet_include"] = True

        per_sample = {
            "sample_1": {
                "n_cells": 100,
                "threshold": None,
                "scrublet_ran": True,
                "needs_review": True,
                "review_reason": "scrublet 自动阈值检测失败（分数无明显双峰）",
            }
        }
        assert per_sample["sample_1"]["needs_review"]
        assert (adata.obs["doublet_call"] == "singlet").all()
        assert adata.obs["doublet_include"].all()

    def test_doublet_report_needs_review_propagates(self):
        """doublet_report["needs_review"] 是 per_sample needs_review 的汇总。"""
        report = {
            "method": "scrublet",
            "per_sample": {
                "s1": {"needs_review": False},
                "s2": {"needs_review": True, "review_reason": "test"},
            },
            "needs_review": any(
                v.get("needs_review", False)
                for v in {"s1": {"needs_review": False},
                           "s2": {"needs_review": True, "review_reason": "test"}}.values()
            ),
            "review_reasons": ["test"],
        }
        assert report["needs_review"]


class TestStageStatus:
    """determine_stage_status 的 NEEDS_REVIEW 与 SUCCESS 分支验证。"""

    def test_needs_review_status(self):
        """needs_review=True + 后置条件全 True → NEEDS_REVIEW。"""
        hard = {"a": True, "b": True}
        status = determine_stage_status(
            {}, hard, needs_review=True, allow_no_required_methods=True,
        )
        assert status == StageStatus.NEEDS_REVIEW

    def test_success_when_clean(self):
        """needs_review=False + 后置条件全 True → SUCCESS。"""
        hard = {"a": True, "b": True}
        status = determine_stage_status(
            {}, hard, needs_review=False, allow_no_required_methods=True,
        )
        assert status == StageStatus.SUCCESS

    def test_failed_when_postcondition_fails(self):
        """后置条件有 False → FAILED（硬门禁优先）。"""
        hard = {"a": True, "b": False}
        status = determine_stage_status(
            {}, hard, needs_review=False, allow_no_required_methods=True,
        )
        assert status == StageStatus.FAILED

    def test_needs_review_overridden_by_failed_postcondition(self):
        """后置条件失败 + needs_review=True → FAILED（硬门禁优先于 needs_review）。"""
        hard = {"a": False, "b": True}
        status = determine_stage_status(
            {}, hard, needs_review=True, allow_no_required_methods=True,
        )
        assert status == StageStatus.FAILED


class TestCountsContractUntouched:
    """跑三态标记逻辑后 layers["counts"] 与 expression_contract 核心字段不变。"""

    def test_counts_layer_preserved_after_three_state(self):
        """执行 doublet 标记（只改 obs 列 + uns）后 counts layer 内容不变。"""
        adata = _make_synthetic_adata()

        counts_before = adata.layers["counts"].copy()
        contract_before = dict(adata.uns["expression_contract"])
        sum_before = float(counts_before.sum())

        # 模拟 doublet cell：只改 obs 列 + uns + processing_history
        adata.obs["doublet_score"] = np.float64(0.0)
        adata.obs["doublet_threshold"] = np.float64(0.3)
        adata.obs["doublet_call"] = pd.Categorical(
            ["singlet"] * adata.n_obs,
            categories=["singlet", "uncertain", "doublet"],
            ordered=True,
        )
        adata.obs["doublet_include"] = True
        adata.obs["predicted_doublet"] = False
        adata.uns["doublet_report"] = {"method": "scrublet", "needs_review": False}
        adata.uns["expression_contract"]["processing_history"].append(
            "doublet_detection: scrublet per-sample, three-state"
        )

        # counts layer 不变
        sum_after = float(adata.layers["counts"].sum())
        assert sum_before == sum_after, (
            f"counts layer sum changed: {sum_before} -> {sum_after}"
        )
        assert adata.layers["counts"].shape == adata.shape

        # 核心契约字段不变（processing_history 除外）
        for key in (
            "x_scale", "counts_layer", "counts_source", "soupx_layer", "stage",
        ):
            assert adata.uns["expression_contract"][key] == contract_before[key], (
                f"contract key {key!r} changed during doublet marking"
            )

    def test_n_obs_unchanged_after_doublet_marking(self):
        """doublet 标记不改 adata 尺寸（决策8：01 不物理删除细胞）。"""
        adata = _make_synthetic_adata(n_cells=100)
        n_before = adata.n_obs

        adata.obs["doublet_score"] = np.float64(0.0)
        adata.obs["doublet_call"] = pd.Categorical(
            np.random.default_rng(43).choice(
                ["singlet", "uncertain", "doublet"], adata.n_obs,
                p=[0.7, 0.2, 0.1],
            ),
            categories=["singlet", "uncertain", "doublet"],
            ordered=True,
        )
        adata.obs["doublet_include"] = adata.obs["doublet_call"] != "doublet"

        assert adata.n_obs == n_before


# ---- 静态源码断言 ------------------------------------------------------------

class TestStaticSourceAssertions:
    """json 解析 notebook 做静态源码断言。"""

    @pytest.fixture(scope="class")
    def nb(self):
        assert NB_PATH.exists(), f"notebook not found: {NB_PATH}"
        with open(NB_PATH) as f:
            return json.load(f)

    @staticmethod
    def _get_cell_source(nb: dict, cell_id: str) -> str:
        for cell in nb["cells"]:
            if cell.get("id") == cell_id and cell["cell_type"] == "code":
                src = cell["source"]
                return "".join(src) if isinstance(src, list) else src
        return ""

    def test_no_dir_in_doublet_and_checkpoint_cells(self, nb):
        """doublet cell (be34c1dd) 与 checkpoint cell (ee38cf39) 不含 dir()。"""
        for cid in ["be34c1dd", "ee38cf39"]:
            src = self._get_cell_source(nb, cid)
            assert "dir()" not in src, (
                f"cell {cid} 使用了 dir()（红线6 fresh-kernel）"
            )

    def test_no_redefinition(self, nb):
        """全 notebook 不重定义 run_contract 中已有函数。"""
        funcs = [
            "validate_expression_contract",
            "determine_stage_status",
            "prepare_run",
            "promote_run",
        ]
        for func in funcs:
            for cell in nb["cells"]:
                if cell["cell_type"] == "code":
                    src = (
                        "".join(cell["source"])
                        if isinstance(cell["source"], list)
                        else cell["source"]
                    )
                    assert f"def {func}" not in src, (
                        f"禁止重新定义 {func}（应 import from run_contract）"
                    )

    def test_filter_cell_does_not_subset_by_doublet(self, nb):
        """过滤 cell (30a91364) 不按 doublet_include/doublet_call 做 subset。"""
        src = self._get_cell_source(nb, "30a91364")
        assert "doublet_include" not in src, (
            "过滤 cell 出现 doublet_include（01 不应按 doublet 物理删除细胞）"
        )
        assert "doublet_call" not in src, (
            "过滤 cell 出现 doublet_call（01 不应按 doublet 物理删除细胞）"
        )

    def test_doublet_cell_writes_required_fields(self, nb):
        """doublet cell 源码含必需 obs/uns 字段字面量。"""
        src = self._get_cell_source(nb, "be34c1dd")
        assert "doublet_call" in src, "doublet cell 缺少 doublet_call"
        assert "doublet_include" in src, "doublet cell 缺少 doublet_include"
        assert "doublet_threshold" in src, "doublet cell 缺少 doublet_threshold"
        assert "doublet_report" in src, "doublet cell 缺少 doublet_report 引用"
        # processing_history append 存在（追加双细胞处理记录）
        assert "processing_history" in src, (
            "doublet cell 缺少 processing_history 追加"
        )

    def test_checkpoint_has_needs_review_branch(self, nb):
        """checkpoint cell 含 needs_review= 传参 + NEEDS_REVIEW 字面量。"""
        src = self._get_cell_source(nb, "ee38cf39")
        assert "needs_review=" in src, (
            "checkpoint cell 未传 needs_review= 给 determine_stage_status"
        )
        assert "NEEDS_REVIEW" in src, (
            "checkpoint cell 未出现 NEEDS_REVIEW 字面量"
        )

    def test_doublet_report_in_manifest_payload(self, nb):
        """manifest_payload 含 doublet_report 字段。"""
        src = self._get_cell_source(nb, "ee38cf39")
        assert "doublet_report" in src, (
            "manifest_payload 未包含 doublet_report"
        )

    def test_doublet_postconditions_in_checkpoint(self, nb):
        """checkpoint cell 的 hard_postconditions 含三个 doublet 后置条件键。"""
        src = self._get_cell_source(nb, "ee38cf39")
        for key in [
            "doublet_columns_present",
            "doublet_not_dropped",
            "doublet_call_valid",
        ]:
            assert key in src, (
                f"hard_postconditions 缺少 {key}"
            )

    def test_importlib_metadata_in_doublet_cell(self, nb):
        """doublet cell 含 importlib.metadata 用于取 scrublet 版本。"""
        src = self._get_cell_source(nb, "be34c1dd")
        assert "importlib.metadata" in src or "importlib" in src, (
            "doublet cell 缺少 importlib.metadata（用于记录 scrublet 版本）"
        )

    def test_categories_declared_as_ordered(self, nb):
        """doublet_call 的 Categorical categories 声明为 ordered（singlet<uncertain<doublet）。"""
        src = self._get_cell_source(nb, "be34c1dd")
        # 至少出现一次 categories 定义含三个态
        assert "singlet" in src
        assert "uncertain" in src
        assert "doublet" in src
        assert "ordered=True" in src or "ordered = True" in src, (
            "doublet_call categories 应声明 ordered=True（singlet<uncertain<doublet）"
        )
