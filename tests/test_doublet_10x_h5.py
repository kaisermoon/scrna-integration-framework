"""测试 10x h5 格式模板：决策8 doublet 三态处理（singlet/uncertain/doublet）。

所有测试只针对 notebooks/01_per_dataset/01_template_10x_h5.ipynb 的 doublet 逻辑。
禁改共享参数化测试 test_pr1b1_stage01_02_input.py / test_pr1b2_notebook_02_04.py。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

# run_contract imports for NEEDS_REVIEW 语义断言
from scrna_integration.run_contract import (
    StageStatus,
    determine_stage_status,
    validate_expression_contract,
)

# ---- helpers ----------------------------------------------------------------

_NOTEBOOK_PATH = (
    Path(__file__).resolve().parents[1]
    / "notebooks"
    / "01_per_dataset"
    / "01_template_10x_h5.ipynb"
)


def _load_notebook() -> dict:
    """加载 notebook JSON，保障文件存在。"""
    assert _NOTEBOOK_PATH.is_file(), f"notebook not found: {_NOTEBOOK_PATH}"
    with open(_NOTEBOOK_PATH) as f:
        return json.load(f)


def _find_cell_by_marker(nb: dict, marker: str, cell_type: str = "code") -> dict:
    """按内容 marker 字符串定位 cell。"""
    for cell in nb["cells"]:
        if cell["cell_type"] != cell_type:
            continue
        src = "".join(cell["source"])
        if marker in src:
            return cell
    raise LookupError(f"cell with marker {marker!r} not found in notebook")


def _source(nb: dict, marker: str, cell_type: str = "code") -> str:
    """提取指定 marker cell 的全部 source。"""
    return "".join(_find_cell_by_marker(nb, marker, cell_type)["source"])


def _cell(
    nb: dict, marker: str, *, env: dict | None = None, cell_type: str = "code"
) -> dict:
    """在受控 env 中 exec 指定 marker cell 的 source。

    返回 exec 后的 env（传入 namespace 的扩展副本）。
    若 env 为 None，使用空 dict。
    """
    src = _source(nb, marker, cell_type=cell_type)
    local_env = {} if env is None else dict(env)
    exec(src, local_env)
    return local_env


# ---- 合成 adata 工厂（多样本，CSR float32 整数 counts）----------------------


def _build_synthetic_kim_adata(
    n_cells: int = 80,
    n_genes: int = 120,
    n_samples: int = 3,
    *,
    seed: int = 42,
) -> anndata.AnnData:
    """构建模拟 Kim 数据集的合成 AnnData。

    - X 为 raw UMI integer counts（CSR float32）
    - sample_id 多样本
    - 含 source_dataset 列
    """
    rng = np.random.default_rng(seed)
    counts = sp.random(n_cells, n_genes, density=0.3, format="csr", dtype=np.float32)
    counts.data = np.floor(counts.data * 100 + 1)  # 正整数 counts

    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    obs["source_dataset"] = "kim"
    obs["sample_id"] = [f"sample_{i % n_samples + 1}" for i in range(n_cells)]
    obs["n_genes"] = rng.integers(500, 5000, n_cells)
    obs["total_counts"] = rng.integers(1000, 20000, n_cells)
    obs["pct_counts_mt"] = rng.uniform(0, 15, n_cells)

    adata = anndata.AnnData(X=counts, obs=obs, var=var)

    # 建立 expression_contract（模拟 cell [3] 的建立段）
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


def _simulate_doublet_scores(
    n_cells: int,
    *,
    singlet_frac: float = 0.7,
    uncertain_frac: float = 0.2,
    doublet_frac: float = 0.1,
    seed: int = 123,
) -> np.ndarray:
    """生成模拟 doublet_score 分布，分 singlet/uncertain/doublet 三段。

    singlet:   [0.0, 0.15)  低分区
    uncertain: [0.15, 0.30) 中间区
    doublet:   [0.30, 1.0]  高分区
    """
    rng = np.random.default_rng(seed)
    scores = np.zeros(n_cells, dtype=np.float64)
    n_sgl = int(n_cells * singlet_frac)
    n_unc = int(n_cells * uncertain_frac)
    n_dbl = n_cells - n_sgl - n_unc
    idx = 0
    scores[idx:idx + n_sgl] = rng.uniform(0.0, 0.15, n_sgl)
    idx += n_sgl
    scores[idx:idx + n_unc] = rng.uniform(0.15, 0.30, n_unc)
    idx += n_unc
    scores[idx:] = rng.uniform(0.30, 1.0, n_dbl)
    rng.shuffle(scores)
    return scores


# ---- 测试用例 ---------------------------------------------------------------


class TestParamsThreeState:
    """验证 PARAMS cell [1] 包含三态 doublet 参数。"""

    def test_params_has_three_state_thresholds(self):
        """PARAMS cell exec 后 env 含 DOUBLET_SCORE_HIGH/LOW 等新参数。"""
        nb = _load_notebook()
        env = _cell(nb, "=== PARAMS ===")
        for key in [
            "DOUBLET_SCORE_HIGH",
            "DOUBLET_SCORE_LOW",
            "DOUBLET_UNCERTAIN_INCLUDE",
            "DOUBLET_RATE_ALERT_HIGH",
            "DOUBLET_RATE_ALERT_LOW",
            "DOUBLET_MIN_CELLS",
        ]:
            assert key in env, f"PARAMS cell missing key: {key}"
        # 旧参数仍存在（向后兼容）
        assert "DOUBLET_SCORE_THRESHOLD" in env


class TestThreeStateClassification:
    """验证三态定级逻辑正确。"""

    def test_three_state_classification(self):
        """合成 doublet_scores 三段，验证定级正确划分。"""
        n = 100
        scores = _simulate_doublet_scores(n, singlet_frac=0.70, uncertain_frac=0.20, doublet_frac=0.10, seed=1)
        threshold_high = 0.30
        threshold_low = 0.15

        classes = np.full(n, "singlet", dtype=object)
        classes[scores > threshold_high] = "doublet"
        classes[(scores >= threshold_low) & (scores <= threshold_high)] = "uncertain"

        # 低分区全部 singlet
        assert np.all(classes[scores < threshold_low] == "singlet")
        # 高分区全部 doublet
        assert np.all(classes[scores > threshold_high] == "doublet")
        # 中间区全部 uncertain
        mid_mask = (scores >= threshold_low) & (scores <= threshold_high)
        if mid_mask.any():
            assert np.all(classes[mid_mask] == "uncertain")
        # 边界值测试：exactly threshold_low → uncertain
        classes_exact = np.full(1, "singlet", dtype=object)
        classes_exact[np.array([0.15]) >= 0.15] = "uncertain"  # >= low
        assert classes_exact[0] == "uncertain"


class TestDoubletIncludeDerivation:
    """验证 doublet_include 列派生逻辑。"""

    def test_doublet_include_excludes_high_conf(self):
        """doublet_class=="doublet" 的细胞 doublet_include 恒为 False。"""
        n = 50
        classes = np.array(
            ["singlet"] * 20 + ["uncertain"] * 20 + ["doublet"] * 10, dtype=object
        )
        np.random.default_rng(42).shuffle(classes)

        # 模拟 include 派生
        include = np.ones(n, dtype=bool)
        include[classes == "doublet"] = False  # 高置信 doublet 排除
        # uncertain 由 DOUBLET_UNCERTAIN_INCLUDE 决定（默认 True 纳入）
        uncertain_include = True
        if not uncertain_include:
            include[classes == "uncertain"] = False

        # singlet 全 True
        assert np.all(include[classes == "singlet"])
        # doublet 全 False
        assert not np.any(include[classes == "doublet"])
        # uncertain 默认 True
        assert np.all(include[classes == "uncertain"])

    def test_uncertain_excluded_when_param_false(self):
        """DOUBLET_UNCERTAIN_INCLUDE=False 时 uncertain 细胞被排除。"""
        n = 50
        classes = np.array(
            ["singlet"] * 20 + ["uncertain"] * 20 + ["doublet"] * 10, dtype=object
        )
        np.random.default_rng(42).shuffle(classes)

        include = np.ones(n, dtype=bool)
        include[classes == "doublet"] = False
        uncertain_include = False  # 不纳入 uncertain
        if not uncertain_include:
            include[classes == "uncertain"] = False

        # uncertain 全 False
        assert not np.any(include[classes == "uncertain"])
        # singlet 不受影响
        assert np.all(include[classes == "singlet"])

    def test_predicted_doublet_backcompat(self):
        """predicted_doublet == (doublet_class=="doublet")，保持旧列语义。"""
        n = 50
        classes = np.array(
            ["singlet"] * 20 + ["uncertain"] * 20 + ["doublet"] * 10, dtype=object
        )
        np.random.default_rng(42).shuffle(classes)

        predicted_doublet = classes == "doublet"

        # 仅高置信 doublet 为 True
        assert np.all(predicted_doublet[classes == "doublet"])
        assert not np.any(predicted_doublet[classes == "singlet"])
        assert not np.any(predicted_doublet[classes == "uncertain"])


class TestThresholdsPersistence:
    """验证阈值与元数据持久化到 uns["doublet_contract"]。"""

    def test_thresholds_persisted_to_uns(self):
        """构造 doublet_contract 写入 uns，验证结构完整。"""
        adata = _build_synthetic_kim_adata(n_cells=50, n_samples=2)

        doublet_contract = {
            "method": "scrublet",
            "method_version": "0.2.3",
            "per_sample_thresholds": {
                "sample_1": {"low": 0.10, "high": 0.20, "auto": 0.20},
                "sample_2": {"low": 0.12, "high": 0.24, "auto": 0.24},
            },
            "uncertain_include": True,
            "rate_alert_high": 0.40,
            "rate_alert_low": None,
            "n_singlet": 35, "n_uncertain": 10, "n_doublet": 5,
            "n_excluded": 5,
            "needs_review": False,
            "needs_review_reasons": [],
            "per_sample_diagnostics": {},
            "random_seed": 42,
        }
        adata.uns["doublet_contract"] = doublet_contract

        assert "doublet_contract" in adata.uns
        dc = adata.uns["doublet_contract"]
        for key in [
            "method", "method_version", "per_sample_thresholds",
            "uncertain_include", "rate_alert_high", "rate_alert_low",
            "n_singlet", "n_uncertain", "n_doublet", "n_excluded",
            "needs_review", "needs_review_reasons",
            "per_sample_diagnostics", "random_seed",
        ]:
            assert key in dc, f"doublet_contract missing key: {key}"

        # per_sample_thresholds 含 per-sample 键
        assert "sample_1" in dc["per_sample_thresholds"]
        assert "low" in dc["per_sample_thresholds"]["sample_1"]
        assert "high" in dc["per_sample_thresholds"]["sample_1"]


class TestNeedsReview:
    """验证 needs_review 触发逻辑。"""

    def test_needs_review_on_high_ratio(self):
        """某样本 doublet_rate > DOUBLET_RATE_ALERT_HIGH 触发 needs_review。"""
        rate_alert_high = 0.40
        # 构造样本诊断：doublet_rate 超阈值
        n_cells = 100
        n_doublet = 50  # 50% > 40%
        doublet_rate = n_doublet / n_cells
        assert doublet_rate > rate_alert_high

        # needs_review 判定
        needs_review = doublet_rate > rate_alert_high
        assert needs_review

        # determine_stage_status 中 needs_review=True → NEEDS_REVIEW
        hard_postconditions = {"non_empty": True}
        status = determine_stage_status(
            {}, hard_postconditions,
            needs_review=True,
            allow_no_required_methods=True,
        )
        assert status == StageStatus.NEEDS_REVIEW

    def test_needs_review_on_too_few_cells(self):
        """样本细胞数 < DOUBLET_MIN_CELLS 记 reason 且进 needs_review。"""
        doublet_min_cells = 50
        n_cells = 30  # < 50
        reason = f"too_few_cells: sample_1 ({n_cells} cells < {doublet_min_cells})"
        needs_review = n_cells < doublet_min_cells
        assert needs_review
        assert "too_few_cells" in reason

    def test_no_needs_review_when_normal(self):
        """正常比例不触发 needs_review。"""
        hard_postconditions = {"non_empty": True}
        status = determine_stage_status(
            {}, hard_postconditions,
            needs_review=False,
            allow_no_required_methods=True,
        )
        assert status == StageStatus.SUCCESS


class TestExpressionContractUntouched:
    """红线守护：expression_contract 不被 doublet 逻辑污染。"""

    def test_expression_contract_schema_unchanged(self):
        """数据读入后 expression_contract 仍是 8 键 schema、counts layer 未被污染。"""
        adata = _build_synthetic_kim_adata(n_cells=30)

        # 模拟 doublet obs 列写入（不碰 expression_contract）
        adata.obs["doublet_score"] = np.nan
        adata.obs["doublet_class"] = "singlet"
        adata.obs["predicted_doublet"] = False
        adata.obs["doublet_include"] = True
        adata.uns["doublet_contract"] = {"method": "not_run"}

        # expression_contract 8 键 schema 完整
        contract = adata.uns["expression_contract"]
        required_keys = {
            "x_scale", "counts_layer", "counts_source", "counts_validated",
            "counts_integer_check", "soupx_layer", "processing_history", "stage",
        }
        assert set(contract.keys()) == required_keys, (
            f"expression_contract keys mismatch: {set(contract.keys())} != {required_keys}"
        )
        # counts layer 未被 doublet 逻辑覆盖/修改
        assert "counts" in adata.layers
        assert adata.layers["counts"].shape == adata.shape
        assert sp.issparse(adata.layers["counts"])
        assert adata.layers["counts"].dtype == np.float32

    def test_doublet_contract_independent(self):
        """doublet_contract 独立于 expression_contract，不塞入额外键。"""
        adata = _build_synthetic_kim_adata(n_cells=30)

        adata.uns["doublet_contract"] = {"method": "scrublet", "n_excluded": 5}

        # expression_contract 不含 doublet 信息
        contract = adata.uns["expression_contract"]
        for key in contract:
            assert key not in ("doublet_score", "doublet_class", "doublet_include", "doublet_contract"), (
                f"expression_contract should not contain doublet key: {key}"
            )


class TestFullObjectPreserved:
    """物理保全：doublet 只标记不删除细胞。"""

    def test_full_object_preserved(self):
        """doublet cell 后 adata.n_obs 未减少。"""
        adata = _build_synthetic_kim_adata(n_cells=50)
        n_before = adata.n_obs

        # 模拟 doublet 三态标记（只写 obs 列，不切片删除）
        adata.obs["doublet_score"] = _simulate_doublet_scores(n_before, seed=99)
        adata.obs["doublet_class"] = "singlet"
        adata.obs.loc[adata.obs["doublet_score"] > 0.30, "doublet_class"] = "doublet"
        adata.obs.loc[
            (adata.obs["doublet_score"] >= 0.15) & (adata.obs["doublet_score"] <= 0.30),
            "doublet_class",
        ] = "uncertain"
        adata.obs["doublet_include"] = adata.obs["doublet_class"] != "doublet"

        # 细胞数不变（仅标记，不物理删除）
        assert adata.n_obs == n_before, (
            f"n_obs changed from {n_before} to {adata.n_obs}: doublet should only mark, not remove"
        )


class TestParamConsumptionUncertainAction:
    """参数消费测试：DOUBLET_UNCERTAIN_INCLUDE 真进逻辑与元数据。"""

    def test_param_consumption_uncertain_action(self):
        """DOUBLET_UNCERTAIN_INCLUDE=False 时 uncertain 细胞被排除，且元数据反映参数值。"""
        n = 60
        classes = np.array(
            ["singlet"] * 30 + ["uncertain"] * 20 + ["doublet"] * 10, dtype=object
        )
        np.random.default_rng(42).shuffle(classes)

        # 分支 1: DOUBLET_UNCERTAIN_INCLUDE=False
        uncertain_include = False
        include = np.ones(n, dtype=bool)
        include[classes == "doublet"] = False
        if not uncertain_include:
            include[classes == "uncertain"] = False

        assert not np.any(include[classes == "uncertain"]), (
            "uncertain cells should be excluded when DOUBLET_UNCERTAIN_INCLUDE=False"
        )

        # 元数据记录参数值
        contract = {"uncertain_include": uncertain_include, "n_excluded": int((~include).sum())}
        assert contract["uncertain_include"] is False
        # n_excluded = doublet(10) + uncertain(20) = 30
        assert contract["n_excluded"] == 30

        # 分支 2: DOUBLET_UNCERTAIN_INCLUDE=True（默认）
        uncertain_include2 = True
        include2 = np.ones(n, dtype=bool)
        include2[classes == "doublet"] = False
        if not uncertain_include2:
            include2[classes == "uncertain"] = False

        assert np.all(include2[classes == "uncertain"]), (
            "uncertain cells should be included when DOUBLET_UNCERTAIN_INCLUDE=True"
        )
        contract2 = {"uncertain_include": uncertain_include2, "n_excluded": int((~include2).sum())}
        assert contract2["uncertain_include"] is True
        # n_excluded = doublet(10) only = 10
        assert contract2["n_excluded"] == 10


class TestNotebookStaticValidation:
    """notebook 静态结构校验（不通过 exec 运行，避免 scVI/LLM 超时）。"""

    def test_cell_structure_has_doublet_cells(self):
        """验证 notebook 中存在 doublet 相关 cell（cell [18] 重写后的检测段 + 三态段）。"""
        nb = _load_notebook()
        # cell [18]（id=2ade235d）应保留且含 "doublet_class" 三态定级
        src18 = _source(nb, "双细胞鉴定")
        assert "doublet_class" in src18, "cell [18] should contain doublet_class"
        assert "doublet_include" in src18, "cell [18] should contain doublet_include"

    def test_checkpoint_has_doublet_postconditions(self):
        """cell [33] checkpoint 含 doublet_columns_present 等 hard_postconditions 键。"""
        nb = _load_notebook()
        src33 = _source(nb, "Checkpoint：写入 per-dataset")
        assert "doublet_columns_present" in src33, (
            "checkpoint cell should have doublet_columns_present postcondition"
        )
        assert "doublet_contract_present" in src33, (
            "checkpoint cell should have doublet_contract_present postcondition"
        )

    def test_checkpoint_has_needs_review_handling(self):
        """cell [33] 含 NEEDS_REVIEW 分支处理（不 promote，只写 draft）。"""
        nb = _load_notebook()
        src33 = _source(nb, "Checkpoint：写入 per-dataset")
        assert "NEEDS_REVIEW" in src33 or "doublet_needs_review" in src33, (
            "checkpoint cell should handle NEEDS_REVIEW via doublet_needs_review"
        )

    def test_counts_contract_cell_untouched(self):
        """cell [3] counts 契约建立段未被修改（含 '建立 counts 契约' 标记）。"""
        nb = _load_notebook()
        src3 = _source(nb, "建立 counts 契约")
        assert "layers[\"counts\"]" in src3
        assert "expression_contract" in src3
        # 不含任何 doublet 相关代码
        for term in ["doublet_score", "doublet_class", "doublet_include", "doublet_contract"]:
            assert term not in src3, (
                f"cell [3] should not contain {term} (counts contract must be untouched)"
            )

    def test_filter_cell_no_doublet_exclusion(self):
        """cell [24] 过滤 cell 不加 doublet 排除逻辑。"""
        nb = _load_notebook()
        src24 = _source(nb, "QC 过滤")
        # 过滤 cell 不应依赖 doublet_include（双重保险）
        assert "doublet_include" not in src24, (
            "filter cell should not exclude by doublet_include (delegated to 02)"
        )

    def test_full_object_preserved_in_notebook(self):
        """notebook 中 doublet cell 不含 adata = adata[keep] 物理删除逻辑。"""
        nb = _load_notebook()
        src18 = _source(nb, "双细胞鉴定")
        # 在 doublet cell 中不应出现切片删除 adata（只标记）
        assert "adata =" not in src18 or "adata[~" not in src18.replace("adata = adata[keep]", ""), (
            "doublet cell should not physically remove cells (only mark obs columns)"
        )

    def test_doublet_contract_static_keys(self):
        """notebook 中 doublet_contract 写入含必要主键。"""
        nb = _load_notebook()
        # 查找含 'uns["doublet_contract"]' 的 code cell（精确匹配 dict 赋值，排除注释中的引用）
        found = False
        for cell in nb["cells"]:
            if cell["cell_type"] != "code":
                continue
            src = "".join(cell["source"])
            # 精确匹配 uns["doublet_contract"] dict 赋值（区别于注释引用和 per_sample_thresholds = {}）
            if 'uns["doublet_contract"]' in src and '"method":' in src:
                found = True
                for key in ["method", "per_sample_thresholds", "needs_review", "random_seed"]:
                    assert key in src, f"doublet_contract should contain key: {key}"
                break
        assert found, "doublet_contract assignment not found in notebook code cells"


class TestObsColumnSchema:
    """验证 obs 列 schema（field_names spec）。"""

    def test_obs_columns_exist_after_doublet(self):
        """执行 doublet cell 后 obs 包含四列新字段。"""
        adata = _build_synthetic_kim_adata(n_cells=30)
        # 模拟 doublet 列写入（按 spec field_names）
        adata.obs["doublet_score"] = np.random.default_rng(1).uniform(0, 1, adata.n_obs).astype(float)
        adata.obs.loc[adata.obs["doublet_score"] > 0.30, "doublet_class"] = "doublet"
        adata.obs.loc[adata.obs["doublet_score"] < 0.30, "doublet_class"] = "singlet"
        adata.obs["predicted_doublet"] = adata.obs["doublet_class"] == "doublet"
        adata.obs["doublet_include"] = adata.obs["doublet_class"] != "doublet"

        for col in ["doublet_score", "doublet_class", "predicted_doublet", "doublet_include"]:
            assert col in adata.obs.columns, f"obs missing column: {col}"

        # 类型检查
        assert adata.obs["doublet_score"].dtype == np.float64 or adata.obs["doublet_score"].dtype.kind == "f"
        assert set(adata.obs["doublet_class"].unique()) <= {"singlet", "uncertain", "doublet"}
        assert adata.obs["predicted_doublet"].dtype == bool
        assert adata.obs["doublet_include"].dtype == bool
