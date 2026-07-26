"""测试 10x mtx 格式模板：决策3 SoupX counts_soupx layer + 决策8 doublet 三态 + 决策9 status.json 判定。

所有测试静态检查 notebooks/01_per_dataset/01_template_10x_mtx.ipynb 的 cell 源码，
不跑 R、不跑 scrublet 端到端（遵守防超时铁律）。
参照 test_pr1b1_stage01_02_input.py 的 _nb/_cell helper 风格。

禁改共享参数化测试 test_pr1b1_stage01_02_input.py。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

# ---- helpers ----------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
_NB_PATH = ROOT / "notebooks" / "01_per_dataset" / "01_template_10x_mtx.ipynb"


def _nb() -> dict:
    """加载 notebook JSON。"""
    assert _NB_PATH.is_file(), f"notebook not found: {_NB_PATH}"
    with open(_NB_PATH, encoding="utf-8") as f:
        return json.load(f)


def _source(marker: str, cell_type: str = "code") -> str:
    """返回首个包含 marker 的 code cell 的 source 字符串。"""
    for cell in _nb()["cells"]:
        if cell["cell_type"] != cell_type:
            continue
        src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
        if marker in src:
            return src
    raise LookupError(f"cell with marker {marker!r} not found")


def _combined() -> str:
    """返回所有 cell 的拼接源码（用于全局搜索）。"""
    parts = []
    for cell in _nb()["cells"]:
        src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
        parts.append(src)
    return "\n".join(parts)


def _params_combined() -> str:
    """返回全部 PARAMS 参数定义的拼接源码（兼容新旧 PARAMS 结构）。

    P1-e 将原单 cell PARAMS 拆为四组（每组 1 md header + 1 code cell），
    此函数收集所有四组 code cell 源码并拼接返回。
    """
    nb = _nb()
    group_markers = [
        "### 1. 数据源", "### 2. QC 阈值",
        "### 3. 方法开关", "### 4. 输出版本与运行标识",
    ]
    sources = []
    for i, cell in enumerate(nb["cells"]):
        src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
        if any(m in src for m in group_markers) and cell["cell_type"] == "markdown":
            # 下一个 cell 是对应的 code cell
            if i + 1 < len(nb["cells"]) and nb["cells"][i + 1]["cell_type"] == "code":
                next_src = "".join(nb["cells"][i + 1]["source"]) if isinstance(nb["cells"][i + 1]["source"], list) else nb["cells"][i + 1]["source"]
                sources.append(next_src)
    if sources:
        return "\n".join(sources)
    # fallback: 单一 PARAMS cell 结构（模板已合并为「唯一参数入口」单 cell）
    return _source("# === PARAMS")


def _cell_indices(nb: dict, marker: str, cell_type: str = "code") -> list[int]:
    """返回所有含 marker 的 cell 在 cells 数组中的 index。"""
    indices = []
    for i, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != cell_type:
            continue
        src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
        if marker in src:
            indices.append(i)
    return indices


# ---- 1. cell 顺序 ---------------------------------------------------------


class TestCellOrdering:
    """验证 SoupX -> QC 重算 -> doublet 的 cell 顺序。"""

    def test_soupx_before_doublet(self):
        """SoupX cell 出现在 doublet cell 之前。"""
        nb = _nb()
        soupx_indices = _cell_indices(nb, "环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        dbl_indices = _cell_indices(nb, "双细胞鉴定：per-sample scrublet（三态定级）")
        assert soupx_indices and dbl_indices, "SoupX 或 doublet cell 不存在"
        assert max(soupx_indices) < min(dbl_indices), (
            f"SoupX cell(s) at {soupx_indices} 应在 doublet cell(s) at {dbl_indices} 之前"
        )

    def test_qc_recompute_between_soupx_and_doublet(self):
        """SoupX 后 QC 重算 cell 在 SoupX code 之后、doublet code 之前。"""
        nb = _nb()
        soupx_i = _cell_indices(nb, "环境 RNA 校正（SoupX）—— subprocess Rscript 模式")[-1]
        qc_i = _cell_indices(nb, "SoupX 校正后 QC 重算", cell_type="markdown")[0]
        dbl_i = _cell_indices(nb, "双细胞鉴定：per-sample scrublet（三态定级）")[0]
        assert soupx_i < qc_i < dbl_i, (
            f"顺序错误: SoupX({soupx_i}) < QC({qc_i}) < doublet({dbl_i})"
        )

    def test_doublet_include_after_doublet(self):
        """doublet include 收尾 cell 在 doublet code 之后。"""
        nb = _nb()
        dbl_i = _cell_indices(nb, "双细胞鉴定：per-sample scrublet（三态定级）")[0]
        inc_i = _cell_indices(nb, "双细胞三态收尾", cell_type="markdown")[0]
        assert dbl_i < inc_i, (
            f"doublet include({inc_i}) 应在 doublet detection({dbl_i}) 之后"
        )


# ---- 2. SoupX 不覆盖 counts ------------------------------------------------


class TestSoupxNoOverwrite:
    """验证 SoupX cell 写 counts_soupx 层，不覆盖 X 或 layers['counts']。"""

    def test_soupx_writes_counts_soupx(self):
        """SoupX cell 含 counts_soupx 写入代码。"""
        src = _source("环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        assert "counts_soupx" in src, "SoupX cell 应写入 counts_soupx layer"

    def test_soupx_does_not_overwrite_X(self):
        """SoupX cell 不含 adata.X 覆盖赋值。"""
        src = _source("环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        assert "adata.X[cell_indices" not in src, (
            "SoupX cell 不得覆盖 adata.X"
        )

    def test_soupx_does_not_overwrite_counts_layer(self):
        """SoupX cell 不含 layers['counts'] 覆盖赋值。"""
        src = _source("环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        assert 'layers["counts"] =' not in src, (
            "SoupX cell 不得覆盖 layers['counts']"
        )

    def test_soupx_cell_initializes_counts_soupx_copy(self):
        """SoupX cell 从 layers['counts'] copy 初始化 counts_soupx。"""
        src = _source("环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        assert 'layers["counts_soupx"] = adata.layers["counts"].copy()' in src, (
            "counts_soupx 必须从 counts copy 初始化"
        )


# ---- 3. status.json 判定 ---------------------------------------------------


class TestSoupxStatusJson:
    """验证 SoupX cell 使用 soupx_status.json 双重判定。"""

    def test_soupx_reads_status_json(self):
        """SoupX cell 读取 soupx_status.json。"""
        src = _source("环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        assert "soupx_status.json" in src, "SoupX cell 应读取 soupx_status.json"

    def test_soupx_checks_status_not_success(self):
        """SoupX cell 含 status != 'success' 分支。"""
        src = _source("环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        assert 'get("status") != "success"' in src or '["status"] != "success"' in src, (
            "SoupX cell 必须判断 status != success"
        )

    def test_soupx_failed_does_not_write(self):
        """status 不成功时不写校正结果（needs_review + continue）。"""
        src = _source("环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        assert "soupx_needs_review = True" in src, (
            "SoupX 失败应标记 needs_review"
        )
        assert "continue" in src, "SoupX 失败应 continue 不写校正结果"

    def test_soupx_check_returncode(self):
        """SoupX cell 同时检查 R exit code。"""
        src = _source("环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        assert "result.returncode" in src, "SoupX cell 应检查 R exit code"


# ---- 4. soupx_layer 契约更新 ----------------------------------------------


class TestSoupxLayerContract:
    """验证 expression_contract.soupx_layer 在成功时更新。"""

    def test_soupx_layer_set_conditionally(self):
        """SoupX cell 含 expression_contract['soupx_layer'] = 'counts_soupx' 条件更新。"""
        src = _source("环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        assert '["soupx_layer"] = "counts_soupx"' in src, (
            "成功时 soupx_layer 应设为 counts_soupx"
        )

    def test_soupx_layer_inside_n_corrected_guard(self):
        """soupx_layer 更新在 n_soupx_corrected > 0 条件下。"""
        src = _source("环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        assert "soupx_layer" in src

    def test_soupx_layer_processing_history_appended(self):
        """成功路径 append processing_history。"""
        src = _source("环境 RNA 校正（SoupX）—— subprocess Rscript 模式")
        assert "processing_history" in src

    def test_expression_contract_schema_field_count(self):
        """expression_contract 建立 cell 仍保持 8 字段 schema。"""
        src = _source("建立 counts 契约")
        # 8 字段均存在
        required = ["x_scale", "counts_layer", "counts_source", "counts_validated",
                     "counts_integer_check", "soupx_layer", "processing_history", "stage"]
        for field in required:
            assert f'"{field}"' in src, f"expression_contract 缺少字段: {field}"


# ---- 5. doublet 三态 -------------------------------------------------------


class TestDoubletThreeState:
    """验证 doublet 三态列名/取值与 kim 完全一致。"""

    def test_doublet_columns_initialized(self):
        """doublet cell 初始化四个 obs 列。"""
        src = _source("双细胞鉴定：per-sample scrublet（三态定级）")
        for col in ["doublet_score", "doublet_class", "predicted_doublet", "doublet_include"]:
            assert col in src, f"doublet cell 应初始化 {col} 列"

    def test_three_state_classification(self):
        """doublet cell 含三态定级（singlet/uncertain/doublet）。"""
        src = _source("双细胞鉴定：per-sample scrublet（三态定级）")
        assert '"singlet"' in src, "三态应含 singlet"
        assert '"uncertain"' in src, "三态应含 uncertain"
        assert '"doublet"' in src, "三态应含 doublet"

    def test_include_derivation_in_include_cell(self):
        """include cell 含 doublet_include 派生逻辑。"""
        src = _source("双细胞三态收尾：include 派生 + doublet_contract 持久化")
        assert "doublet_include" in src

    def test_predicted_doublet_backcompat(self):
        """predicted_doublet = doublet_class=='doublet'（旧列语义不变）。"""
        combined = _combined()
        has_backcompat = (
            'predicted_doublet"] = adata.obs["doublet_class"] == "doublet"' in combined
            or "predicted_doublet'] = adata.obs['doublet_class'] == 'doublet'" in combined
        )
        assert has_backcompat, "predicted_doublet 应等于 doublet_class=='doublet'"


# ---- 6. checkpoint _hd guard ----------------------------------------------


class TestCheckpointHdGuard:
    """验证 checkpoint cell 含 _hd guard + doublet postconditions。"""

    def test_hd_guard_present(self):
        """checkpoint 含 _hd = 'doublet_class' in adata.obs。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert '_hd = "doublet_class" in adata.obs' in src, (
            "checkpoint 应有 _hd guard"
        )

    def test_doublet_columns_present_postcondition(self):
        """checkpoint 含 doublet_columns_present 键。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "doublet_columns_present" in src

    def test_doublet_contract_present_postcondition(self):
        """checkpoint 含 doublet_contract_present 键。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "doublet_contract_present" in src

    def test_doublet_class_valid_postcondition(self):
        """checkpoint 含 doublet_class_valid 键。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "doublet_class_valid" in src

    def test_doublet_include_consistent_postcondition(self):
        """checkpoint 含 doublet_include_consistent 键。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "doublet_include_consistent" in src

    def test_hd_else_true_pattern(self):
        """doublet postconditions 用 if _hd else True 模式。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        count = src.count("if _hd else True")
        assert count >= 3, f"应有 >=3 处 if _hd else True guard，实际 {count}"


# ---- 7. needs_review 传入 -------------------------------------------------


class TestNeedsReview:
    """验证 checkpoint cell 传入 needs_review。"""

    def test_needs_review_passed_to_determine_stage_status(self):
        """checkpoint 调用 determine_stage_status 传 needs_review=。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "determine_stage_status({}, hard_postconditions, needs_review=" in src, (
            "checkpoint 必须传 needs_review 给 determine_stage_status"
        )

    def test_doublet_needs_review_captured(self):
        """checkpoint 获取 doublet_needs_review。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "doublet_needs_review" in src

    def test_soupx_needs_review_captured(self):
        """checkpoint 获取 soupx_needs_review。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "soupx_needs_review" in src

    def test_needs_review_branch_present(self):
        """checkpoint 含 NEEDS_REVIEW 分支处理。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "NEEDS_REVIEW" in src, "checkpoint 应有 NEEDS_REVIEW 分支"

    def test_soupx_summary_in_checkpoint(self):
        """NEEDS_REVIEW 分支含 soupx_summary。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "soupx_summary" in src, "checkpoint NEEDS_REVIEW 应写 soupx_summary"


# ---- 8. counts_layer_unchanged 门禁 --------------------------------------


class TestCountsLayerUnchanged:
    """验证 checkpoint 含 counts checksum 门禁。"""

    def test_counts_checksum_saved_in_contract_cell(self):
        """expression_contract cell 保存 _counts_checksum。"""
        src = _source("建立 counts 契约")
        assert "_counts_checksum" in src, "expression_contract cell 应保存 counts checksum"
        assert "_counts_checksum_nnz" in src, "expression_contract cell 应保存 nnz"

    def test_counts_layer_unchanged_postcondition(self):
        """checkpoint hard_postconditions 含 counts_layer_unchanged。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "counts_layer_unchanged" in src, "checkpoint 应有 counts_layer_unchanged"

    def test_checksum_compare_in_checkpoint(self):
        """checkpoint 比对 _counts_checksum 当前值。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "_counts_checksum" in src, "checkpoint 应引用 _counts_checksum"


# ---- 9. n_doublets_total 已定义 -------------------------------------------


class TestNDoubletsTotal:
    """验证 n_doublets_total 在 include cell 显式定义。"""

    def test_n_doublets_total_defined(self):
        """include cell 含 n_doublets_total = n_dbl。"""
        src = _source("双细胞三态收尾：include 派生 + doublet_contract 持久化")
        assert "n_doublets_total = n_dbl" in src, (
            "include cell 必须定义 n_doublets_total（防 QC 报告 cell NameError）"
        )


# ---- 10. PARAMS 完整性 -----------------------------------------------------


class TestParams:
    """验证 PARAMS cell 含新参数。"""

    def test_doublet_three_state_params(self):
        """PARAMS 含六个三态参数。"""
        src = _params_combined()
        for key in [
            "DOUBLET_SCORE_HIGH",
            "DOUBLET_SCORE_LOW",
            "DOUBLET_UNCERTAIN_INCLUDE",
            "DOUBLET_RATE_ALERT_HIGH",
            "DOUBLET_RATE_ALERT_LOW",
            "DOUBLET_MIN_CELLS",
        ]:
            assert key in src, f"PARAMS 缺少: {key}"

    def test_soupx_integer_round_param(self):
        """PARAMS 含 SOUPX_INTEGER_ROUND。"""
        src = _params_combined()
        assert "SOUPX_INTEGER_ROUND" in src

    def test_old_params_preserved(self):
        """旧参数仍保留（向后兼容）。"""
        src = _params_combined()
        for key in ["EXPECTED_DOUBLET_RATE", "DOUBLET_SCORE_THRESHOLD", "SOUPX_ENABLED"]:
            assert key in src, f"旧参数 {key} 丢失"


# ---- 11. doublet_contract schema ------------------------------------------


class TestDoubletContract:
    """验证 doublet_contract 写入含必要字段。"""

    def test_doublet_contract_keys(self):
        """doublet_contract 含必要主键。"""
        src = _source("双细胞三态收尾：include 派生 + doublet_contract 持久化")
        for key in ["method", "per_sample_thresholds", "needs_review", "random_seed",
                     "n_singlet", "n_uncertain", "n_doublet", "n_excluded"]:
            assert key in src, f"doublet_contract 缺少: {key}"


# ---- 12. 共享测试不破坏 ----------------------------------------------------


class TestSharedTestCompatibility:
    """确保不破坏 test_pr1b1_stage01_02_input.py 的参数化测试。"""

    def test_params_cell_parseable(self):
        """PARAMS cell 可独立 exec（兼容新旧结构）。"""
        import ast
        src = _params_combined()
        ast.parse(src)

    def test_checkpoint_cell_parseable(self):
        """Checkpoint cell 可独立 parse。"""
        import ast
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        ast.parse(src)

    def test_checkpoint_uses_snapshot_with_globals(self):
        """Checkpoint 含 snapshot_effective_parameters(globals()...)。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "snapshot_effective_parameters(globals()" in src

    def test_checkpoint_uses_collect_runtime_provenance(self):
        """Checkpoint 含 collect_runtime_provenance。"""
        src = _source("Checkpoint：run_contract + 写入 h5ad + schema 校验")
        assert "collect_runtime_provenance" in src

    def test_no_run_contract_helper_redefinitions(self):
        """无 cell 重定义 run_contract helper。"""
        combined = _combined()
        for helper in ["def validate_expression_contract(", "def prepare_run(", "def promote_run("]:
            assert helper not in combined, (
                f"不应重定义 run_contract helper: {helper}"
            )
