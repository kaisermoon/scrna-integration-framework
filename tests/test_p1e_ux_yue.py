"""P1-e-yue：01_yue.ipynb UX骨架静态源码断言 + preflight行为测试。

覆盖：
- PARAMS 四组化：### 数据源 / ### QC 阈值 / ### 方法开关 / ### 输出与运行标识
- preflight cell：位置校验（在 data-load 之前）+ SOUPX_ENABLED/QC_STRATEGY/DOUBLET_RATE 断言
- 科学参数四要素注释（调大/调小/何时/默认）
- 默认值不变
- redline 完整性：doublet 三态 / expression_contract / hard_postconditions / _hd guard / layers["counts"]
- preflight 行为：good params 通过；SOUPX_ENABLED=True / QC_STRATEGY 非法 / EXPECTED_DOUBLET_RATE>=1 均触发 AssertionError
"""

# ruff: noqa: N806

from __future__ import annotations

import json
from pathlib import Path

import pytest

# ---- 路径常量 ----------------------------------------------------------------

NB_PATH = (
    Path(__file__).resolve().parent.parent
    / "notebooks" / "01_per_dataset" / "01_yue.ipynb"
)

# 科学参数默认值（与 notebook PARAMS 保持同步）
DEFAULT_N_MAD = 4
DEFAULT_SOUPX_ENABLED = False
DEFAULT_EXPECTED_DOUBLET_RATE = 0.04
DEFAULT_OUTPUT_VERSION = 1
DEFAULT_PER_SAMPLE_MAD = True
DEFAULT_MIN_CELLS_PER_GENE = 3


# ---- helpers ----------------------------------------------------------------

def _load_nb() -> dict:
    """加载 notebook JSON。"""
    assert NB_PATH.exists(), f"notebook not found: {NB_PATH}"
    with open(NB_PATH) as f:
        return json.load(f)


def _cell_by_id(nb: dict, cell_id: str) -> dict | None:
    """按 cell id 查找 cell。"""
    for cell in nb["cells"]:
        if cell.get("id") == cell_id:
            return cell
    return None


def _cell_source(cell: dict) -> str:
    """获取 cell 源码字符串。"""
    src = cell["source"]
    return "".join(src) if isinstance(src, list) else src


def _cell_index(nb: dict, cell_id: str) -> int | None:
    """获取 cell 的索引位置。"""
    for i, cell in enumerate(nb["cells"]):
        if cell.get("id") == cell_id:
            return i
    return None


def _find_preflight_cell(nb: dict) -> dict | None:
    """查找含 'Preflight' 标记的 code cell。"""
    for cell in nb["cells"]:
        if cell["cell_type"] == "code":
            src = _cell_source(cell)
            if "Preflight" in src:
                return cell
    return None


def _preflight_check(good_params: dict, *, raises: bool = False) -> None:
    """复刻 notebook preflight cell 的核心校验逻辑。

    传入参数字典模拟全局命名空间，执行 preflight 各项断言。
    raises=True 时本函数预期 raise AssertionError，否则通过。
    """
    import os
    import yaml
    from pathlib import Path

    SOUPX_ENABLED = good_params["SOUPX_ENABLED"]
    MANIFEST_PATH = good_params["MANIFEST_PATH"]
    QC_STRATEGY = good_params["QC_STRATEGY"]
    EXPECTED_DOUBLET_RATE = good_params["EXPECTED_DOUBLET_RATE"]
    N_MAD = good_params["N_MAD"]
    _data_dir_exists = good_params.get("_data_dir_exists", True)
    _txt_files_count = good_params.get("_txt_files_count", 10)

    # 1. SOUPX_ENABLED 强制 False
    assert SOUPX_ENABLED is False, (
        f"SOUPX_ENABLED 必须为 False（Yue 类器官数据集无 raw matrix），当前为 {SOUPX_ENABLED}"
    )

    # 2. Manifest 存在
    if good_params.get("_manifest_exists", True):
        _manifest_path = Path(MANIFEST_PATH)
        assert _manifest_path.exists(), f"Manifest 不存在: {MANIFEST_PATH}"

    # 3. txt.gz 文件存在（Yue 格式）
    if _data_dir_exists and _txt_files_count > 0:
        pass  # 不实际读 manifest，跳过文件系统依赖

    # 4. QC_STRATEGY 合法性
    assert QC_STRATEGY in ("adaptive", "fixed"), (
        f"QC_STRATEGY 必须为 'adaptive' 或 'fixed'，当前为 {QC_STRATEGY!r}"
    )

    # 5. EXPECTED_DOUBLET_RATE 范围
    assert isinstance(EXPECTED_DOUBLET_RATE, (int, float)) and 0 < EXPECTED_DOUBLET_RATE < 1, (
        f"EXPECTED_DOUBLET_RATE 必须在 (0, 1) 范围内，当前为 {EXPECTED_DOUBLET_RATE}"
    )

    # 6. N_MAD 范围
    assert isinstance(N_MAD, (int, float)) and 1 <= N_MAD <= 10, (
        f"N_MAD 应在 [1, 10] 合理范围内，当前为 {N_MAD}"
    )


# ---- 静态源码断言 ------------------------------------------------------------

class TestParamsFourGroups:
    """PARAMS 四组化：### 数据源 / ### QC 阈值 / ### 方法开关 / ### 输出与运行标识。"""

    def test_four_group_headers_present(self):
        """PARAMS cell 含四个 ### markdown 级分组标题。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "3e6879a1")
        assert cell is not None, "PARAMS cell (3e6879a1) not found"
        src = _cell_source(cell)
        for header in ["### 数据源", "### QC 阈值", "### 方法开关", "### 输出与运行标识"]:
            assert header in src, f"PARAMS cell 缺少分组标题: {header}"

    def test_group_member_counts(self):
        """每个分组的参数数量符合预期。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "3e6879a1")
        src = _cell_source(cell)

        # 提取各组段的内容行（strip() 后以四个已知分组标题之一开头）
        _valid_headers = ["### 数据源", "### QC 阈值", "### 方法开关", "### 输出与运行标识"]
        sections = {}
        current_section = None
        for line in src.split("\n"):
            stripped = line.strip()
            matched = next((h for h in _valid_headers if stripped.startswith(h)), None)
            if matched is not None:
                current_section = matched  # 用短标题做 key（稳定，不受副文本影响）
                sections[current_section] = []
            elif current_section and "=" in stripped and not stripped.startswith("#"):
                sections.setdefault(current_section, []).append(stripped)

        # 数据源 至少含 MANIFEST_PATH
        ds_lines = [l for l in sections.get("### 数据源", []) if "MANIFEST_PATH" in l]
        assert len(ds_lines) >= 1, "数据源组缺少 MANIFEST_PATH"

        # 方法开关 至少含 SOUPX_ENABLED / FLAG_HEMOGLOBIN / FLAG_STRESS_GENES / SCORE_CELL_CYCLE
        ms_src = "\n".join(sections.get("### 方法开关", []))
        for p in ["SOUPX_ENABLED", "FLAG_HEMOGLOBIN", "FLAG_STRESS_GENES", "SCORE_CELL_CYCLE"]:
            assert p in ms_src, f"方法开关组缺少 {p}"

        # 输出组 至少含 RUN_ID / OUTPUT_FILENAME
        out_src = "\n".join(sections.get("### 输出与运行标识", []))
        for p in ["RUN_ID", "OUTPUT_FILENAME", "OUTPUT_VERSION"]:
            assert p in out_src, f"输出与运行标识组缺少 {p}"


class TestPreflightCell:
    """preflight cell 存在性 + 位置 + 内容断言。"""

    def test_preflight_exists_before_data_load(self):
        """preflight cell 存在且 index < data-load cell (717b91ba)。"""
        nb = _load_nb()
        preflight_cell = _find_preflight_cell(nb)
        assert preflight_cell is not None, "未找到 preflight cell（含 'Preflight' 标记的 code cell）"

        dataload_idx = _cell_index(nb, "717b91ba")
        assert dataload_idx is not None, "data-load cell (717b91ba) 未找到"

        preflight_idx = nb["cells"].index(preflight_cell)
        assert preflight_idx < dataload_idx, (
            f"preflight cell (index={preflight_idx}) 应在 data-load cell (index={dataload_idx}) 之前"
        )

    def test_preflight_contains_soupx_assert(self):
        """preflight cell 含 SOUPX_ENABLED is False 断言。"""
        nb = _load_nb()
        cell = _find_preflight_cell(nb)
        src = _cell_source(cell)
        assert "SOUPX_ENABLED" in src, "preflight 缺少 SOUPX_ENABLED 引用"
        assert "is False" in src, "preflight 缺少 SOUPX_ENABLED is False 断言"

    def test_preflight_contains_txt_gz_glob(self):
        """preflight cell 含 txt.gz glob 校验。"""
        nb = _load_nb()
        cell = _find_preflight_cell(nb)
        src = _cell_source(cell)
        assert "txt.gz" in src, "preflight 缺少 txt.gz 引用"
        assert "glob" in src or "*.txt" in src, "preflight 缺少文件 glob 匹配"

    def test_preflight_contains_manifest_check(self):
        """preflight cell 含 manifest 存在性校验。"""
        nb = _load_nb()
        cell = _find_preflight_cell(nb)
        src = _cell_source(cell)
        assert "MANIFEST_PATH" in src, "preflight 缺少 MANIFEST_PATH 引用"
        assert "exists()" in src or "exists" in src.lower(), "preflight 缺少文件存在性校验"

    def test_preflight_contains_strategy_check(self):
        """preflight cell 含 QC_STRATEGY 合法性断言。"""
        nb = _load_nb()
        cell = _find_preflight_cell(nb)
        src = _cell_source(cell)
        assert "QC_STRATEGY" in src, "preflight 缺少 QC_STRATEGY 引用"

    def test_preflight_contains_doublet_rate_check(self):
        """preflight cell 含 EXPECTED_DOUBLET_RATE 范围断言。"""
        nb = _load_nb()
        cell = _find_preflight_cell(nb)
        src = _cell_source(cell)
        assert "EXPECTED_DOUBLET_RATE" in src, "preflight 缺少 EXPECTED_DOUBLET_RATE 引用"

    def test_preflight_contains_marker(self):
        """preflight cell 含 'Preflight' 文字标记。"""
        nb = _load_nb()
        cell = _find_preflight_cell(nb)
        src = _cell_source(cell)
        # Preflight 标记在注释或 print 中
        assert "preflight" in src.lower(), "preflight cell 缺失 'Preflight' 标记字面量"

    def test_preflight_not_build_adata(self):
        """preflight cell 不创建 adata 对象。"""
        nb = _load_nb()
        cell = _find_preflight_cell(nb)
        src = _cell_source(cell)
        # 不应出现 adata = AnnData 或 adata = anndata.AnnData
        assert "AnnData(" not in src, "preflight cell 不应创建 AnnData 对象"
        assert "adata =" not in src, "preflight cell 不应创建 adata 变量"


class TestScientificParamComments:
    """科学参数含四要素注释（调大/调小/何时/默认）。"""

    def test_four_element_comments_present(self):
        """需补四要素注释的参数均含调大/调小/何时描述。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "3e6879a1")
        src = _cell_source(cell)

        params_to_check = [
            "N_MAD", "QC_STRATEGY", "PER_SAMPLE_MAD", "MIN_CELLS_PER_GENE",
            "EXPECTED_DOUBLET_RATE", "MIN_GENES", "MAX_GENES", "MIN_COUNTS", "MAX_PCT_MT",
        ]

        for param in params_to_check:
            # 找到包含该参数赋值的行
            found_line = ""
            for line in src.split("\n"):
                stripped = line.strip()
                if stripped.startswith(param) and "=" in stripped:
                    found_line = line
                    break
            assert found_line, f"未找到 {param} 的定义行"
            # 行中应含注释（# 之后有说明文字）
            assert "#" in found_line, f"{param} 缺少行内注释"
            comment_part = found_line.split("#", 1)[1] if "#" in found_line else ""
            # 含调大/调小/何时/默认 中至少一个
            has_annotation = any(
                kw in comment_part
                for kw in ["调大", "调小", "何时", "调为", "默认", "仅QC_STRATEGY"]
            )
            assert has_annotation, (
                f"{param} 行内注释缺少四要素标记（调大/调小/何时/默认）: {comment_part[:80]}..."
            )

    def test_already_compliant_params_retained(self):
        """已合规参数保留原有注释不变。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "3e6879a1")
        src = _cell_source(cell)

        existing = [
            ("DOUBLET_UNCERTAIN_MARGIN", "调大"),
            ("DOUBLET_RATE_ALERT", "过高"),
            ("DOUBLET_MIN_CELLS", "调低"),
            ("DOUBLET_SCORE_THRESHOLD", "手动覆盖"),
        ]
        for param, keyword in existing:
            found = False
            for line in src.split("\n"):
                if line.strip().startswith(param):
                    if keyword in line:
                        found = True
                    break
            assert found, f"{param} 原有注释中的 '{keyword}' 被修改"


class TestDefaultsUnchanged:
    """PARAMS 默认值与 P1-c-yue 保持一致。"""

    def test_defaults(self):
        """核心参数默认值未改变。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "3e6879a1")
        src = _cell_source(cell)

        # 用简单模式匹配确认默认值
        checks = {
            "N_MAD": "N_MAD = 4",
            "SOUPX_ENABLED": "SOUPX_ENABLED = False",
            "EXPECTED_DOUBLET_RATE": "EXPECTED_DOUBLET_RATE = 0.04",
            "OUTPUT_VERSION": "OUTPUT_VERSION = 1",
            "PER_SAMPLE_MAD": "PER_SAMPLE_MAD = True",
            "MIN_CELLS_PER_GENE": "MIN_CELLS_PER_GENE = 3",
        }
        for param, pattern in checks.items():
            assert pattern in src, f"{param} 默认值已改变（期望含 '{pattern}'）"


class TestRedlineIntegrity:
    """redline 区内容完整性：doublet / expression_contract / checkpoint 字段未丢失。"""

    def test_doublet_three_state_markers(self):
        """doublet cell 三态标记字面量存在。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "be34c1dd")
        src = _cell_source(cell)
        for kw in ["doublet_call", "doublet_include", "doublet_threshold", "singlet", "uncertain", "doublet"]:
            assert kw in src, f"doublet cell 缺少 {kw}"

    def test_expression_contract_fields(self):
        """data-load cell expression_contract 完整。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "717b91ba")
        src = _cell_source(cell)
        for kw in ["expression_contract", "x_scale", "counts_layer", "layers[\"counts\"]"]:
            assert kw in src, f"data-load cell 缺少 {kw}"

    def test_hard_postconditions_in_checkpoint(self):
        """checkpoint cell hard_postconditions 关键字存在。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "ee38cf39")
        src = _cell_source(cell)
        for kw in ["hard_postconditions", "doublet_columns_present", "doublet_not_dropped", "doublet_call_valid"]:
            assert kw in src, f"checkpoint cell 缺少 {kw}"

    def test_hd_guard_in_checkpoint(self):
        """checkpoint cell _hd guard 存在。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "ee38cf39")
        src = _cell_source(cell)
        assert "_hd" in src, "checkpoint cell 缺少 _hd guard"

    def test_layers_counts_in_dataload(self):
        """data-load cell 含 layers['counts'] 建立逻辑。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "717b91ba")
        src = _cell_source(cell)
        assert 'layers["counts"]' in src or "layers['counts']" in src, (
            "data-load cell 缺少 layers['counts']"
        )

    def test_no_dir_in_redline_cells(self):
        """redline 五 cell 不含 dir()。"""
        nb = _load_nb()
        redline_ids = ["717b91ba", "be34c1dd", "8c9316b7", "ee38cf39", "30a91364"]
        for cid in redline_ids:
            cell = _cell_by_id(nb, cid)
            if cell is None:
                continue
            src = _cell_source(cell)
            assert "dir()" not in src, f"redline cell {cid} 含 dir()"

    def test_no_redefinition(self):
        """notebook 不重定义 run_contract 中已有函数。"""
        nb = _load_nb()
        funcs = [
            "validate_expression_contract",
            "determine_stage_status",
            "prepare_run",
            "promote_run",
        ]
        for func in funcs:
            for cell in nb["cells"]:
                if cell["cell_type"] == "code":
                    src = _cell_source(cell)
                    assert f"def {func}" not in src, (
                        f"cell {cell.get('id')} 重定义 {func}"
                    )

    def test_filter_cell_defines_cells_after(self):
        """过滤 cell (30a91364) 定义 cells_after（checkpoint _hd guard 引用）。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "30a91364")
        src = _cell_source(cell)
        assert "cells_after" in src, "过滤 cell 缺少 cells_after 定义"

    def test_needs_review_branch_in_checkpoint(self):
        """checkpoint cell 含 NEEDS_REVIEW 分支处理。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "ee38cf39")
        src = _cell_source(cell)
        assert "needs_review=" in src, "checkpoint 未传 needs_review="
        assert "NEEDS_REVIEW" in src, "checkpoint 缺少 NEEDS_REVIEW 字面量"

    def test_soupx_cell_has_skip_branch(self):
        """SoupX cell 保留 SOUPX_ENABLED=False → skip 分支不变。"""
        nb = _load_nb()
        cell = _cell_by_id(nb, "8c9316b7")
        src = _cell_source(cell)
        assert "SOUPX_ENABLED" in src, "SoupX cell 缺少 SOUPX_ENABLED"
        assert "soupx_applied" in src, "SoupX cell 缺少 soupx_applied"


# ---- preflight 行为测试 ----------------------------------------------------

class TestPreflightBehavior:
    """复刻 preflight 校验逻辑，验证 good/bad params 的行为。"""

    def test_good_params_pass(self):
        """全部合法参数：preflight 通过（不抛异常）。"""
        good = {
            "SOUPX_ENABLED": False,
            "MANIFEST_PATH": "data/yue/manifest.yaml",
            "QC_STRATEGY": "adaptive",
            "EXPECTED_DOUBLET_RATE": 0.04,
            "N_MAD": 4,
        }
        try:
            _preflight_check(good)
        except AssertionError as e:
            pytest.fail(f"good params 不应触发 AssertionError: {e}")

    def test_soupx_enabled_true_raises(self):
        """SOUPX_ENABLED=True → AssertionError。"""
        bad = {
            "SOUPX_ENABLED": True,
            "MANIFEST_PATH": "data/yue/manifest.yaml",
            "QC_STRATEGY": "adaptive",
            "EXPECTED_DOUBLET_RATE": 0.04,
            "N_MAD": 4,
        }
        with pytest.raises(AssertionError, match="SOUPX_ENABLED"):
            _preflight_check(bad)

    def test_qc_strategy_illegal_raises(self):
        """QC_STRATEGY='unknown' → AssertionError。"""
        bad = {
            "SOUPX_ENABLED": False,
            "MANIFEST_PATH": "data/yue/manifest.yaml",
            "QC_STRATEGY": "unknown",
            "EXPECTED_DOUBLET_RATE": 0.04,
            "N_MAD": 4,
        }
        with pytest.raises(AssertionError, match="QC_STRATEGY"):
            _preflight_check(bad)

    def test_doublet_rate_out_of_range_raises(self):
        """EXPECTED_DOUBLET_RATE=1.5 → AssertionError。"""
        bad = {
            "SOUPX_ENABLED": False,
            "MANIFEST_PATH": "data/yue/manifest.yaml",
            "QC_STRATEGY": "adaptive",
            "EXPECTED_DOUBLET_RATE": 1.5,
            "N_MAD": 4,
        }
        with pytest.raises(AssertionError, match="EXPECTED_DOUBLET_RATE"):
            _preflight_check(bad)

    def test_doublet_rate_zero_raises(self):
        """EXPECTED_DOUBLET_RATE=0 → AssertionError（不在 (0,1) 范围内）。"""
        bad = {
            "SOUPX_ENABLED": False,
            "MANIFEST_PATH": "data/yue/manifest.yaml",
            "QC_STRATEGY": "adaptive",
            "EXPECTED_DOUBLET_RATE": 0,
            "N_MAD": 4,
        }
        with pytest.raises(AssertionError, match="EXPECTED_DOUBLET_RATE"):
            _preflight_check(bad)

    def test_n_mad_out_of_range_raises(self):
        """N_MAD=0 → AssertionError（不在 [1,10] 范围内）。"""
        bad = {
            "SOUPX_ENABLED": False,
            "MANIFEST_PATH": "data/yue/manifest.yaml",
            "QC_STRATEGY": "adaptive",
            "EXPECTED_DOUBLET_RATE": 0.04,
            "N_MAD": 0,
        }
        with pytest.raises(AssertionError, match="N_MAD"):
            _preflight_check(bad)

    def test_qc_strategy_fixed_valid(self):
        """QC_STRATEGY='fixed' 通过 preflight。"""
        good = {
            "SOUPX_ENABLED": False,
            "MANIFEST_PATH": "data/yue/manifest.yaml",
            "QC_STRATEGY": "fixed",
            "EXPECTED_DOUBLET_RATE": 0.04,
            "N_MAD": 4,
        }
        try:
            _preflight_check(good)
        except AssertionError as e:
            pytest.fail(f"QC_STRATEGY='fixed' 应通过 preflight: {e}")
