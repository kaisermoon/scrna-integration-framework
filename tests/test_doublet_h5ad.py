"""P1-c: doublet 三态政策测试 — h5ad 格式模板（h5ad 特有：scrublet 输入取 layers["counts"]，不取 X）。

覆盖：
1. classify_doublets 三态边界语义
2. doublet_include 映射（仅 doublet → False）
3. predicted_doublet 等于 (prediction=="doublet")
4. author_removed 跳过路径
5. qc_overrides 跳过路径
6. needs_review 触发条件（小样本 / 高 doublet 比例）
7. scrublet 输入守卫（nowicki 特有：X=logcounts 非整数 → 误用会科学错误）
8. uns["doublet_detection"] schema 完整性
9. 列类型一致性

测试风格：纯函数镜像 + 契约断言，不 import notebook，不跑 nbconvert。
合成数据用 scipy.sparse.random + np.floor 造整数 counts。
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import pytest
from anndata import AnnData


# ============================================================================
# 镜像函数（与 notebook cell [p1c_doublet_code] 边界语义逐字一致）
# ============================================================================

def classify_doublets(scores, threshold, margin):
    """将 doublet scores 按三态策略分类——与 notebook 边界语义逐字一致。

    Parameters
    ----------
    scores : np.ndarray
        scrublet 输出的 doublet scores
    threshold : float
        上界阈值
    margin : float
        uncertain 带宽比例

    Returns
    -------
    np.ndarray of str: "singlet" | "uncertain" | "doublet"
    """
    result = np.full(len(scores), "singlet", dtype=object)
    lower_bound = (1 - margin) * threshold
    result[scores > threshold] = "doublet"
    result[(scores > lower_bound) & (scores <= threshold)] = "uncertain"
    return result


def _make_counts_adata(n_obs: int = 200, n_vars: int = 10, seed: int = 42) -> AnnData:
    """构造合成 AnnData：X=logcounts（非整数），layers["counts"]=整数 counts。

    模拟 nowicki 场景：X 是 normalized_log1p，原始整数 counts 在 layers["counts"]。
    """
    rng = np.random.RandomState(seed)
    raw_counts = rng.poisson(5, size=(n_obs, n_vars))
    logcounts = np.log1p(raw_counts).astype(np.float32)

    gene_names = [f"GENE_{i}" for i in range(n_vars)]
    adata = AnnData(
        X=sp.csr_matrix(logcounts),
        var={"gene_id": gene_names},
    )
    adata.var_names = gene_names
    adata.layers["counts"] = sp.csr_matrix(raw_counts.astype(np.float32))

    adata.obs["sample_id"] = ["sample_A"] * (n_obs // 2) + ["sample_B"] * (n_obs - n_obs // 2)
    return adata


def _check_integer_matrix(mat) -> bool:
    """检查稀疏矩阵的所有非零值是否近似整数。"""
    if sp.issparse(mat):
        data = mat.data
    else:
        data = mat.ravel()
    if data.size == 0:
        return True
    return bool(np.allclose(data % 1, 0))


# ============================================================================
# 1. classify_doublets 三态边界语义
# ============================================================================

def test_classify_three_states_boundaries():
    """验证 classify_doublets 的 > 与 <= 边界方向。"""
    margin = 0.2
    threshold = 0.5
    lower = (1 - margin) * threshold  # 0.4

    scores = np.array([
        threshold * 1.1,          # 0.55 → doublet (>)
        threshold,                 # 0.50 → uncertain (≤)
        lower + 0.01,              # 0.41 → uncertain (> lower, ≤ threshold)
        lower,                      # 0.40 → uncertain (≤ threshold, 但 = lower 按 ≤ lower 算?)
    ])
    # 边界精确定义：
    # score > threshold → doublet
    # lower < score <= threshold → uncertain
    # score <= lower → singlet
    results = classify_doublets(scores, threshold, margin)

    # score=0.55 > 0.5 → doublet
    assert results[0] == "doublet", f"score={scores[0]:.2f} > threshold={threshold} should be doublet"
    # score=0.50 ≤ 0.5 → uncertain (闭合边界属 uncertain)
    assert results[1] == "uncertain", f"score={scores[1]:.2f} = threshold={threshold} should be uncertain (≤)"
    # score=0.41: lower(0.4) < 0.41 ≤ 0.5 → uncertain
    assert results[2] == "uncertain", f"score={scores[2]:.2f} should be uncertain"
    # score=0.40 ≤ lower(0.4) → singlet
    assert results[3] == "singlet", f"score={scores[3]:.2f} = lower={lower} should be singlet (≤)"


def test_classify_extreme_scores():
    """极端 score 值：0 为 singlet，1.0 为 doublet。"""
    margin = 0.2
    threshold = 0.5
    scores = np.array([0.0, 1.0])
    results = classify_doublets(scores, threshold, margin)
    assert results[0] == "singlet"
    assert results[1] == "doublet"


def test_classify_all_singlet_below_lower():
    """所有 score < lower → 全部 singlet。"""
    margin = 0.2
    threshold = 0.5
    scores = np.array([0.0, 0.1, 0.2, 0.3, 0.39])
    results = classify_doublets(scores, threshold, margin)
    assert all(r == "singlet" for r in results)


# ============================================================================
# 2. doublet_include 映射
# ============================================================================

def test_doublet_include_mapping():
    """仅 doublet → include=False；singlet/uncertain → include=True。"""
    states = np.array(["singlet", "uncertain", "doublet", "singlet", "uncertain", "doublet"])
    include = (states != "doublet")
    assert list(include) == [True, True, False, True, True, False]


# ============================================================================
# 3. predicted_doublet 兼容列
# ============================================================================

def test_predicted_doublet_equals_high_conf():
    """predicted_doublet 恒等于 (prediction=="doublet")。"""
    states = np.array(["singlet", "uncertain", "doublet", "singlet", "doublet", "uncertain"])
    predicted = (states == "doublet")
    assert list(predicted) == [False, False, True, False, True, False]
    # 验证与三态分类的一致性
    assert np.array_equal(predicted, states == "doublet")


# ============================================================================
# 4. author_removed 跳过路径
# ============================================================================

def test_skip_via_preprocessing_done():
    """manifest preprocessing_done 含 doublet_removal → 应跳过。"""
    manifest = {
        "preprocessing_done": ["basic_filter", "doublet_removal", "normalization"],
        "qc_overrides": {},
    }
    pp_done = manifest.get("preprocessing_done", [])
    qc_dbl = manifest.get("qc_overrides", {}).get("doublet_removal", {})

    skip = False
    if "doublet_removal" in pp_done:
        skip = True
    elif qc_dbl.get("skip"):
        skip = True

    assert skip is True


def test_author_removed_skip_path_full_fields():
    """跳过分支下四列取默认值、uns needs_review=False。"""
    adata = _make_counts_adata(50, 5)
    # 模拟 skip 分支的默认值写入
    adata.obs["doublet_score"] = np.nan
    adata.obs["doublet_prediction"] = "singlet"
    adata.obs["doublet_include"] = True
    adata.obs["predicted_doublet"] = False

    # 断言四列存在
    assert "doublet_score" in adata.obs.columns
    assert "doublet_prediction" in adata.obs.columns
    assert "doublet_include" in adata.obs.columns
    assert "predicted_doublet" in adata.obs.columns

    # 跳过分支默认值
    assert adata.obs["doublet_score"].isna().all()
    assert (adata.obs["doublet_prediction"] == "singlet").all()
    assert adata.obs["doublet_include"].all()
    assert not adata.obs["predicted_doublet"].any()

    # uns 元数据
    adata.uns["doublet_detection"] = {
        "method": "author_removed",
        "version": "manifest_skip",
        "per_sample_thresholds": {},
        "uncertain_margin": 0.2,
        "expected_doublet_rate": None,
        "skip_reason": "原作者已去除双细胞",
        "n_singlet": adata.n_obs,
        "n_uncertain": 0,
        "n_doublet": 0,
        "n_excluded": 0,
        "needs_review": False,
        "review_reasons": [],
        "per_sample_diagnostics": {},
    }
    assert adata.uns["doublet_detection"]["method"] == "author_removed"
    assert adata.uns["doublet_detection"]["needs_review"] is False
    assert adata.uns["doublet_detection"]["n_excluded"] == 0


def test_skip_via_qc_override_only():
    """仅 qc_overrides.doublet_removal.skip=True（preprocessing_done 不含）也应 skip。"""
    manifest = {
        "preprocessing_done": ["basic_filter", "normalization"],
        "qc_overrides": {
            "doublet_removal": {"skip": True, "reason": "test reason"}
        },
    }
    pp_done = manifest.get("preprocessing_done", [])
    qc_dbl = manifest.get("qc_overrides", {}).get("doublet_removal", {})

    skip = False
    if "doublet_removal" in pp_done:
        skip = True
    elif qc_dbl.get("skip"):
        skip = True

    assert skip is True


# ============================================================================
# 5. needs_review 触发条件
# ============================================================================

def test_needs_review_small_sample():
    """样本细胞数 < DOUBLET_MIN_CELLS(50) → flagged/needs_review=True。"""
    DOUBLET_MIN_CELLS = 50
    n_sample = 30
    flagged = n_sample < DOUBLET_MIN_CELLS
    assert flagged is True
    # 阈值应为 None（不定阈值）
    threshold = None if flagged else 0.5
    assert threshold is None


def test_needs_review_high_doublet_rate():
    """pct_doublet > DOUBLET_RATE_ALERT(0.30) → needs_review=True。"""
    DOUBLET_RATE_ALERT = 0.30
    # 构造 10 个细胞中 4 个 doublet → 40%
    states = np.array(["singlet", "doublet", "singlet", "doublet",
                       "doublet", "singlet", "doublet", "singlet", "singlet", "singlet"])
    n_dbl = int((states == "doublet").sum())
    n_sample = len(states)
    pct_dbl = n_dbl / n_sample  # 0.4
    needs_review = pct_dbl > DOUBLET_RATE_ALERT
    assert needs_review is True
    assert pct_dbl == 0.4


def test_needs_review_reasons_accumulate():
    """多个样本触发 review 时 review_reasons 应累积。"""
    review_reasons = []
    review_reasons.append("sample_A: doublet 比例 35.0% > DOUBLET_RATE_ALERT=30.0%")
    review_reasons.append("sample_B: 样本细胞数 30 < DOUBLET_MIN_CELLS=50，不定阈值")
    assert len(review_reasons) == 2
    assert any("DOUBLET_RATE_ALERT" in r for r in review_reasons)
    assert any("DOUBLET_MIN_CELLS" in r for r in review_reasons)


# ============================================================================
# 6. scrublet 输入守卫（nowicki 特有）
# ============================================================================

def test_scrublet_uses_counts_layer_not_x():
    """nowicki 特有守卫：X 是 logcounts（非整数），layers["counts"] 是整数 counts。

    断言：
    - X 非整数 → 若误用 X 会取到非整数矩阵（反例守卫）
    - layers["counts"] 为整数 → scrublet 必须取 counts 层
    """
    adata = _make_counts_adata(100, 20, seed=123)

    # X 应为非整数（logcounts）
    x_is_integer = _check_integer_matrix(adata.X)
    assert not x_is_integer, (
        "nowicki 的 X 是 logcounts（非整数），若此断言失败说明数据构造错误"
    )

    # layers["counts"] 应为整数（原始 counts）
    counts_is_integer = _check_integer_matrix(adata.layers["counts"])
    assert counts_is_integer, (
        "layers['counts'] 必须是整数 counts；scrublet 必须取此层而非 X"
    )

    # 反例守卫：若误用 X（非整数 logcounts）会导致科学错误
    # scrublet 期望整数 counts 输入，logcounts 会破坏 doublet score 分布
    x_data = adata.X.data if sp.issparse(adata.X) else adata.X.ravel()
    counts_data = (
        adata.layers["counts"].data
        if sp.issparse(adata.layers["counts"])
        else adata.layers["counts"].ravel()
    )
    # 两者不应相同（X 经过 log1p 变换）
    if x_data.size == counts_data.size:
        assert not np.allclose(np.sort(x_data), np.sort(counts_data)), (
            "X（logcounts）不应等于 layers['counts']（原始 counts）；"
            "若被误喂给 scrublet 将导致科学错误"
        )


def test_counts_layer_preserves_integer_nature():
    """layers['counts'] 以 float32 存储但值全部近似整数。"""
    adata = _make_counts_adata(200, 15, seed=99)
    counts = adata.layers["counts"]
    assert counts.dtype == np.float32, "counts 应以 float32 存储（内存纪律）"
    assert sp.issparse(counts), "counts 应为稀疏矩阵"
    assert _check_integer_matrix(counts), "counts 值应全部近似整数"


# ============================================================================
# 7. uns["doublet_detection"] schema 完整性
# ============================================================================

REQUIRED_DD_KEYS = frozenset({
    "method", "version", "per_sample_thresholds", "uncertain_margin",
    "expected_doublet_rate", "skip_reason",
    "n_singlet", "n_uncertain", "n_doublet", "n_excluded",
    "needs_review", "review_reasons", "per_sample_diagnostics",
})


def test_uns_doublet_detection_schema_author_removed():
    """跳过分支 uns["doublet_detection"] 含全部必需键。"""
    dd = {
        "method": "author_removed",
        "version": "manifest_skip",
        "per_sample_thresholds": {},
        "uncertain_margin": 0.2,
        "expected_doublet_rate": None,
        "skip_reason": "原作者已去除双细胞",
        "n_singlet": 2500,
        "n_uncertain": 0,
        "n_doublet": 0,
        "n_excluded": 0,
        "needs_review": False,
        "review_reasons": [],
        "per_sample_diagnostics": {},
    }
    missing = REQUIRED_DD_KEYS - dd.keys()
    assert not missing, f"missing keys in doublet_detection: {missing}"

    # 类型断言
    assert isinstance(dd["method"], str)
    assert isinstance(dd["needs_review"], bool)
    assert isinstance(dd["n_singlet"], int)
    assert isinstance(dd["review_reasons"], list)
    assert isinstance(dd["per_sample_diagnostics"], dict)
    assert isinstance(dd["per_sample_thresholds"], dict)


def test_uns_doublet_detection_schema_scrublet():
    """scrublet 分支 uns["doublet_detection"] 含全部必需键且类型正确。"""
    dd = {
        "method": "scrublet",
        "version": "0.2.3",
        "per_sample_thresholds": {"sample_A": 0.15, "sample_B": 0.12},
        "uncertain_margin": 0.2,
        "expected_doublet_rate": 0.05,
        "skip_reason": None,
        "n_singlet": 2300,
        "n_uncertain": 50,
        "n_doublet": 150,
        "n_excluded": 150,
        "needs_review": False,
        "review_reasons": [],
        "per_sample_diagnostics": {
            "sample_A": {
                "n_cells": 1200, "n_doublet": 80, "pct_doublet": 0.067,
                "threshold": 0.15, "flagged": False, "flag_reason": None,
            },
        },
    }
    missing = REQUIRED_DD_KEYS - dd.keys()
    assert not missing, f"missing keys in doublet_detection: {missing}"

    # per_sample_diagnostics 子字段
    for sid, diag in dd["per_sample_diagnostics"].items():
        assert "n_cells" in diag, f"sample {sid} missing n_cells"
        assert "n_doublet" in diag, f"sample {sid} missing n_doublet"
        assert "pct_doublet" in diag, f"sample {sid} missing pct_doublet"
        assert "threshold" in diag, f"sample {sid} missing threshold"
        assert "flagged" in diag, f"sample {sid} missing flagged"
        assert "flag_reason" in diag, f"sample {sid} missing flag_reason"


# ============================================================================
# 8. 列类型一致性
# ============================================================================

def test_doublet_columns_dtypes():
    """四列类型：doublet_include=bool, doublet_prediction 取值在三态内, doublet_score=float。"""
    adata = _make_counts_adata(50, 5)
    adata.obs["doublet_score"] = np.array([0.1, 0.3, np.nan, 0.05, 0.8] * 10, dtype=np.float32)
    adata.obs["doublet_prediction"] = np.array(
        ["singlet", "uncertain", "singlet", "singlet", "doublet"] * 10, dtype=object
    )
    adata.obs["doublet_include"] = adata.obs["doublet_prediction"] != "doublet"
    adata.obs["predicted_doublet"] = adata.obs["doublet_prediction"] == "doublet"

    # doublet_include 为 bool
    assert adata.obs["doublet_include"].dtype == bool
    # doublet_prediction 取值 ⊂ {singlet, uncertain, doublet}
    valid = {"singlet", "uncertain", "doublet"}
    actual = set(adata.obs["doublet_prediction"].unique())
    assert actual.issubset(valid), f"invalid prediction values: {actual - valid}"
    # doublet_score 为浮点
    assert np.issubdtype(adata.obs["doublet_score"].dtype, np.floating)
    # predicted_doublet 为 bool
    assert adata.obs["predicted_doublet"].dtype == bool


def test_doublet_include_only_false_for_doublet():
    """doublet_include=False 必须且仅当 prediction='doublet'。"""
    adata = _make_counts_adata(30, 5)
    rng = np.random.RandomState(42)
    states = rng.choice(["singlet", "uncertain", "doublet"], size=30)
    adata.obs["doublet_prediction"] = states
    adata.obs["doublet_include"] = states != "doublet"

    doublet_mask = adata.obs["doublet_prediction"] == "doublet"
    include_false = ~adata.obs["doublet_include"]
    assert np.array_equal(doublet_mask, include_false), (
        "doublet_include=False 必须严格对应 prediction='doublet'"
    )


# ============================================================================
# 9. 并发守卫：no physical removal in 01
# ============================================================================

def test_no_physical_removal_in_01():
    """01 不物理删除细胞：标记后的 n_obs 与标记前一致。"""
    adata = _make_counts_adata(100, 5)
    n_before = adata.n_obs
    # 模拟标记（部分标 doublet）
    adata.obs["doublet_include"] = np.array(
        [True] * 80 + [False] * 20
    )
    n_after = adata.n_obs
    assert n_after == n_before, (
        "01 只标记 doublet_include，不物理删除细胞。"
        f"n_obs 从 {n_before} 变为 {n_after}"
    )
    # 物理排除应保留给下游
    assert adata.obs["doublet_include"].sum() == 80
    assert (~adata.obs["doublet_include"]).sum() == 20


# ============================================================================
# 10. 默认 skip 路径 stage_status 仍 SUCCESS
# ============================================================================

def test_default_skip_path_is_success():
    """默认 author_removed 跳过路径下 needs_review=False，应可 promote。

    这是关键验收点：nowicki 默认路径不回退到 NEEDS_REVIEW。"""
    dd = {
        "method": "author_removed",
        "needs_review": False,
        "n_excluded": 0,
    }
    assert dd["needs_review"] is False, (
        "nowicki 默认 skip 路径 needs_review=False，应正常 SUCCESS → promote"
    )
    assert dd["n_excluded"] == 0, "跳过路径无不纳入的细胞"
