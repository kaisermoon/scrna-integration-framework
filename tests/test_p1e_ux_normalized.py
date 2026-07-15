"""P1-e-normalized: 03_normalized.ipynb PARAMS 四组化 + preflight 结构合规测试。

覆盖：
- 8 个 PARAMS 四组 cell 存在且 ID 正确、顺序正确
- 2 个 preflight cell 存在且 ID 正确、位于 03-setup 与 03-load 之间
- 所有 code cell AST 可解析
- PARAMS 四组中各变量的右值与原 03-params 逐字一致（只重组、不改语义）
- 受保护 cell（03-load / 03-counts-layer / 03-normalize-code / 03-checkpoint-code）
  源码与 main 原版逐字一致
- preflight cell 含四段校验关键词（参数类型 / manifest / expression_contract / 生效参数）
- 全 notebook 无 dir() 结果判断
- 不在 notebook 中重定义 run_contract 函数
"""

from __future__ import annotations

import ast
import json
import re
import subprocess
from pathlib import Path

import pytest


# ---- 路径常量 -------------------------------------------------------------------

_NB_PATH = Path(__file__).resolve().parent.parent / "notebooks" / "03_normalized.ipynb"

# 原 03-params cell 的变量名（用于值不变性比对）
_ALL_PARAM_NAMES = [
    "UPSTREAM_RUN_ROOT", "UPSTREAM_RUN_ID",
    "NORMALIZATION_METHOD", "TARGET_SUM",
    "N_TOP_GENES", "HVG_FLAVOR",
    "BATCH_AWARE_HVG", "HVG_BATCH_KEY",
    "EXCLUDE_MT_FROM_HVG", "EXCLUDE_RIBO_FROM_HVG", "EXCLUDE_HB_FROM_HVG",
    "CUSTOM_EXCLUDE_PATTERNS",
    "HVG_SUBSAMPLE_PER_BATCH",
    "EXCLUDE_CELL_CYCLE_FROM_HVG",
    "FORCED_INCLUDE_GENES",
    "REGRESS_OUT",
    "SCALE", "MAX_SCALE_VALUE",
    "RUN_ID", "RUN_ROOT", "OUTPUT_FILENAME", "OUTPUT_VERSION",
    "RANDOM_SEED",
]

# 受保护 cell（源码不得有任何改动）
_PROTECTED_CELL_IDS = [
    "03-load", "03-counts-layer", "03-normalize-code", "03-checkpoint-code", "03-setup",
]

# run_contract 函数（notebook 禁止重定义）
_RUN_CONTRACT_FUNCTIONS = [
    "validate_expression_contract", "aggregate_method_status",
    "determine_stage_status", "prepare_run", "promote_run", "resume_run",
    "atomic_write_json", "sha256_file", "validate_checkpoint",
    "collect_runtime_provenance", "snapshot_effective_parameters",
]


# ---- 辅助函数 -------------------------------------------------------------------


def _load_nb() -> dict:
    """加载当前 notebook JSON。"""
    with open(_NB_PATH, encoding="utf-8") as f:
        return json.load(f)


def _load_orig_nb() -> dict:
    """从 git 加载 main 的原始 notebook JSON。"""
    result = subprocess.run(
        ["git", "show", f"HEAD:notebooks/03_normalized.ipynb"],
        capture_output=True, text=True, check=True,
        cwd=str(_NB_PATH.parent.parent),
    )
    return json.loads(result.stdout)


def _get_cell(nb: dict, cell_id: str) -> dict | None:
    """按 id 查找 cell。"""
    for cell in nb["cells"]:
        if cell.get("id") == cell_id:
            return cell
    return None


def _get_cell_source(nb: dict, cell_id: str) -> str:
    """获取 cell 的源码字符串。"""
    cell = _get_cell(nb, cell_id)
    if cell is None:
        return ""
    src = cell.get("source", "")
    return "".join(src) if isinstance(src, list) else src


def _extract_param_values(nb: dict) -> dict[str, str]:
    """从 PARAMS 四组 cell 中提取变量名到右值字符串的映射。"""
    values: dict[str, str] = {}
    param_cell_ids = [
        "03-params-data", "03-params-science", "03-params-compute", "03-params-output",
    ]
    for cid in param_cell_ids:
        src = _get_cell_source(nb, cid)
        for line in src.split("\n"):
            line_stripped = line.strip()
            if line_stripped.startswith("#") or not line_stripped:
                continue
            m = re.match(r"^(\w+)\s*=\s*(.+)$", line_stripped)
            if m:
                name = m.group(1)
                val = m.group(2).strip()
                # 剥去行内注释（保留原值字符串）
                val = re.sub(r"#.*$", "", val).rstrip()
                values[name] = val
    return values


def _all_code_cells(nb: dict):
    """生成所有 code cell 的 (id, source_string)。"""
    for cell in nb["cells"]:
        if cell.get("cell_type") != "code":
            continue
        src = cell.get("source", "")
        yield cell.get("id", "<no-id>"), "".join(src) if isinstance(src, list) else src


# ---- 结构合规测试 ---------------------------------------------------------------


class TestParamsRestructure:
    """PARAMS 四组化结构：8 个新 cell 的存在、ID 与顺序。"""

    _EXPECTED_ORDER = [
        "03-title",
        "03-params-md", "03-params-data",
        "03-params-science-md", "03-params-science",
        "03-params-compute-md", "03-params-compute",
        "03-params-output-md", "03-params-output",
        "03-setup",
    ]

    def test_old_params_cell_removed(self) -> None:
        """原 03-params cell 已删除。"""
        nb = _load_nb()
        assert _get_cell(nb, "03-params") is None, "03-params cell 应已删除"

    def test_eight_new_params_cells_exist(self) -> None:
        """8 个新 PARAMS cell 全部存在且 ID 正确。"""
        nb = _load_nb()
        expected_ids = [
            "03-params-md", "03-params-data",
            "03-params-science-md", "03-params-science",
            "03-params-compute-md", "03-params-compute",
            "03-params-output-md", "03-params-output",
        ]
        for cid in expected_ids:
            cell = _get_cell(nb, cid)
            assert cell is not None, f"缺少 cell: {cid}"

    def test_params_cells_in_correct_order(self) -> None:
        """PARAMS cell 间的相对顺序正确。"""
        nb = _load_nb()
        cell_ids = [c.get("id", "") for c in nb["cells"]]

        # 03-title 和 03-setup 之间存在 8 个新 cell
        title_idx = cell_ids.index("03-title")
        setup_idx = cell_ids.index("03-setup")

        between = cell_ids[title_idx + 1 : setup_idx]
        expected_md_code_alternating = [
            "03-params-md", "03-params-data",
            "03-params-science-md", "03-params-science",
            "03-params-compute-md", "03-params-compute",
            "03-params-output-md", "03-params-output",
        ]
        assert between == expected_md_code_alternating, (
            f"PARAMS cell 顺序不匹配:\n  实际: {between}\n  期望: {expected_md_code_alternating}"
        )

    def test_params_md_cells_are_markdown(self) -> None:
        """PARAMS 组的 markdown cell 确实是 markdown 类型。"""
        nb = _load_nb()
        md_ids = [
            "03-params-md", "03-params-science-md",
            "03-params-compute-md", "03-params-output-md",
        ]
        for cid in md_ids:
            cell = _get_cell(nb, cid)
            assert cell is not None
            assert cell["cell_type"] == "markdown", f"{cid} 应为 markdown，实际为 {cell['cell_type']}"

    def test_params_code_cells_are_code(self) -> None:
        """PARAMS 组的 code cell 确实是 code 类型。"""
        nb = _load_nb()
        code_ids = [
            "03-params-data", "03-params-science",
            "03-params-compute", "03-params-output",
        ]
        for cid in code_ids:
            cell = _get_cell(nb, cid)
            assert cell is not None
            assert cell["cell_type"] == "code", f"{cid} 应为 code，实际为 {cell['cell_type']}"


class TestPreflightStructure:
    """Preflight cell 的存在、ID、位置与关键词。"""

    def test_preflight_cells_exist(self) -> None:
        """两个 preflight cell 存在。"""
        nb = _load_nb()
        assert _get_cell(nb, "03-preflight-md") is not None
        assert _get_cell(nb, "03-preflight-code") is not None

    def test_preflight_between_setup_and_load(self) -> None:
        """preflight cell 位于 03-setup 与 03-load 之间。"""
        nb = _load_nb()
        cell_ids = [c.get("id", "") for c in nb["cells"]]
        setup_idx = cell_ids.index("03-setup")
        load_idx = cell_ids.index("03-load")
        assert cell_ids[setup_idx + 1] == "03-preflight-md", "03-setup 后应为 03-preflight-md"
        assert cell_ids[setup_idx + 2] == "03-preflight-code", "03-preflight-md 后应为 03-preflight-code"
        assert cell_ids[setup_idx + 3] == "03-load", "03-preflight-code 后应为 03-load"

    def test_preflight_md_is_markdown(self) -> None:
        """03-preflight-md 是 markdown cell。"""
        nb = _load_nb()
        cell = _get_cell(nb, "03-preflight-md")
        assert cell["cell_type"] == "markdown"

    def test_preflight_code_is_code(self) -> None:
        """03-preflight-code 是 code cell。"""
        nb = _load_nb()
        cell = _get_cell(nb, "03-preflight-code")
        assert cell["cell_type"] == "code"

    def test_preflight_has_param_type_check(self) -> None:
        """preflight code 含参数类型校验（步骤1）。"""
        src = _get_cell_source(_load_nb(), "03-preflight-code")
        assert "NORMALIZATION_METHOD" in src
        assert "_valid_methods" in src
        assert "TARGET_SUM" in src
        assert "N_TOP_GENES" in src
        assert "HVG_FLAVOR" in src
        assert "_valid_flavors" in src
        assert "isinstance" in src, "应使用 isinstance 做类型检查"

    def test_preflight_has_manifest_check(self) -> None:
        """preflight code 含 manifest 存在性与合法性校验（步骤2）。"""
        src = _get_cell_source(_load_nb(), "03-preflight-code")
        assert "manifest.json" in src
        assert "UPSTREAM_RUN_ID" in src
        assert "02_merged" in src, "应校验 manifest.stage == '02_merged'"

    def test_preflight_has_expression_contract_check(self) -> None:
        """preflight code 含 backed 模式 expression_contract 校验（步骤3）。"""
        src = _get_cell_source(_load_nb(), "03-preflight-code")
        assert 'backed="r"' in src or "backed='r'" in src, "应使用 backed='r' 模式只读元数据"
        assert "expression_contract" in src
        assert "validate_expression_contract" in src
        assert "raw_counts" in src, "应校验 expected_scale='raw_counts'"
        # 校验 layers["counts"] 键存在（backed 下只查键名）
        # preflight 中写法为 if "counts" not in _handle.layers: raise KeyError(...)
        assert '"counts"' in src, "应校验 layers['counts'] 键存在"
        assert "_handle.layers" in src, "应通过 backed handle 检查 layers"
        assert ".file.close()" in src, "backed 模式必须显式关闭 _handle.file.close()"

    def test_preflight_has_echo(self) -> None:
        """preflight code 含生效参数回显（步骤4）。"""
        src = _get_cell_source(_load_nb(), "03-preflight-code")
        assert "生效参数" in src
        assert "NORMALIZATION_METHOD" in src
        assert "RUN_ID" in src

    def test_preflight_has_pass_message(self) -> None:
        """preflight code 末尾打印通过消息。"""
        src = _get_cell_source(_load_nb(), "03-preflight-code")
        assert "Preflight 通过" in src

    def test_preflight_no_soupx_doublet_empty_checks(self) -> None:
        """preflight 不为 SoupX/doublet 加空校验（那是 01 的职责）。
        允许注释提及 doublet/SoupX 说明 N/A 原因，但不允许实际校验代码。
        """
        src = _get_cell_source(_load_nb(), "03-preflight-code")
        # 不应有 soupx_layer 字段检查（那是 01 的契约字段，03 不检查）
        assert "soupx_layer" not in src
        # 允许注释说明 N/A，但不应检查 doublet_prediction / doublet_include 等列
        assert "doublet_prediction" not in src
        assert "doublet_include" not in src
        assert "doublet_class" not in src


# ---- 值不变性测试 ----------------------------------------------------------------


class TestParamsValuePreservation:
    """PARAMS 四组化后所有变量右值与原始 03-params 逐字一致。"""

    def test_all_param_values_identical(self) -> None:
        """每个变量在新四组 cell 中的值与原 03-params cell 完全一致。"""
        orig_nb = _load_orig_nb()
        new_nb = _load_nb()

        # 从原 03-params cell 提取值
        orig_src = ""
        for cell in orig_nb["cells"]:
            if cell.get("id") == "03-params":
                orig_src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
                break
        assert orig_src, "原始 notebook 中未找到 id=03-params"

        orig_values: dict[str, str] = {}
        for line in orig_src.split("\n"):
            line_stripped = line.strip()
            if line_stripped.startswith("#") or not line_stripped:
                continue
            m = re.match(r"^(\w+)\s*=\s*(.+)$", line_stripped)
            if m:
                name = m.group(1)
                val = m.group(2).strip()
                val = re.sub(r"#.*$", "", val).rstrip()
                orig_values[name] = val

        # 从新四组 cell 提取值
        new_values = _extract_param_values(new_nb)

        mismatches = []
        for name in sorted(_ALL_PARAM_NAMES):
            orig_val = orig_values.get(name)
            new_val = new_values.get(name)
            if orig_val != new_val:
                mismatches.append(
                    f"  {name}: orig={orig_val!r} new={new_val!r}"
                )

        # 无缺失
        for name in _ALL_PARAM_NAMES:
            if name not in new_values:
                mismatches.append(f"  {name}: 在新 PARAMS cell 中缺失")

        assert not mismatches, (
            f"PARAMS 值不一致 ({len(mismatches)} 项):\n" + "\n".join(mismatches)
        )

    def test_param_count_match(self) -> None:
        """新四组 cell 包含的变量数与原始一致。"""
        new_nb = _load_nb()
        new_values = _extract_param_values(new_nb)
        assert len(new_values) == len(_ALL_PARAM_NAMES), (
            f"变量数不匹配：新={len(new_values)} vs 原={len(_ALL_PARAM_NAMES)}"
        )


# ---- 受保护 cell 不变性测试 -------------------------------------------------------


class TestProtectedCellsUntouched:
    """受保护 cell 的源码与 main 逐字一致。"""

    def test_protected_cells_unchanged(self) -> None:
        """所有受保护 cell 源码与原始逐字一致。"""
        orig_nb = _load_orig_nb()
        new_nb = _load_nb()

        violations = []
        for cid in _PROTECTED_CELL_IDS:
            orig_cell = _get_cell(orig_nb, cid)
            new_cell = _get_cell(new_nb, cid)
            if orig_cell is None:
                violations.append(f"  {cid}: 原始 notebook 中缺失")
                continue
            if new_cell is None:
                violations.append(f"  {cid}: 新 notebook 中缺失")
                continue
            orig_src = "".join(orig_cell["source"]) if isinstance(orig_cell["source"], list) else orig_cell["source"]
            new_src = "".join(new_cell["source"]) if isinstance(new_cell["source"], list) else new_cell["source"]
            if orig_src != new_src:
                violations.append(f"  {cid}: 源码已变更 ({len(orig_src)} -> {len(new_src)} bytes)")

        assert not violations, (
            f"受保护 cell 不应被修改 ({len(violations)} 项):\n" + "\n".join(violations)
        )


# ---- notebook 静态合规测试 --------------------------------------------------------


class TestNotebookStaticCompliance:
    """通用静态合规检查。"""

    def test_all_code_cells_ast_parse(self) -> None:
        """所有 code cell AST 可解析。"""
        nb = _load_nb()
        errors = []
        for cid, src in _all_code_cells(nb):
            try:
                ast.parse(src)
            except SyntaxError as e:
                errors.append(f"  {cid}: {e}")
        assert not errors, f"AST 解析失败 ({len(errors)}):\n" + "\n".join(errors)

    def test_no_dir_in_any_code_cell(self) -> None:
        """所有 code cell 不使用 dir() 判断变量存在性（红线6）。"""
        nb = _load_nb()
        violations = []
        for cid, src in _all_code_cells(nb):
            if "dir()" in src:
                violations.append(f"  {cid}")
        assert not violations, (
            f"以下 cell 使用了 dir()（红线6）:\n" + "\n".join(violations)
        )

    def test_no_redefinition_of_run_contract_functions(self) -> None:
        """禁止在 notebook 中重定义 run_contract 函数（红线5）。"""
        nb = _load_nb()
        violations = []
        for cid, src in _all_code_cells(nb):
            for func_name in _RUN_CONTRACT_FUNCTIONS:
                if f"def {func_name}" in src:
                    violations.append(f"  {cid}: def {func_name}")
        # 允许 setup cell 中的 import 行
        actual_violations = [v for v in violations if "03-setup" not in v]
        assert not actual_violations, (
            f"以下 cell 重定义了 run_contract 函数:\n" + "\n".join(actual_violations)
        )

    def test_imports_validate_expression_contract_in_setup(self) -> None:
        """03-setup cell 导入了 validate_expression_contract。"""
        src = _get_cell_source(_load_nb(), "03-setup")
        assert "validate_expression_contract" in src, (
            "03-setup 缺少 validate_expression_contract 导入"
        )

    def test_checkpoint_has_counts_layer_preserved(self) -> None:
        """03-checkpoint-code cell 的 hard_postconditions 含 counts_layer_preserved。"""
        src = _get_cell_source(_load_nb(), "03-checkpoint-code")
        assert '"counts_layer_preserved"' in src
        assert "_counts_integrity_checked" in src

    def test_protected_cells_no_dir_usage(self) -> None:
        """受保护 cell 不使用 dir()（红线6 针对性检查）。"""
        for cid in ["03-counts-layer", "03-normalize-code", "03-checkpoint-code"]:
            src = _get_cell_source(_load_nb(), cid)
            assert "dir()" not in src, f"{cid} 使用了 dir()"

    def test_json_legal(self) -> None:
        """notebook JSON 可直接 json.load 解析。"""
        nb = _load_nb()
        assert isinstance(nb, dict)
        assert "cells" in nb
        assert len(nb["cells"]) > 0


# ---- 科学参数四要素注释测试 -------------------------------------------------------


class TestScienceParamAnnotations:
    """03-params-science cell 中每个科学参数上方有四要素注释。"""

    def test_four_element_annotations_present(self) -> None:
        """科学参数 code cell 含四要素注释关键词。"""
        src = _get_cell_source(_load_nb(), "03-params-science")
        # 每个科学参数上方应有：是什么 / 默认依据 / 调大调小 / 何时该改
        assert "是什么" in src, "缺少四要素注释「是什么」"
        assert "默认依据" in src, "缺少四要素注释「默认依据」"
        assert "何时该改" in src, "缺少四要素注释「何时该改」"

    def test_major_params_have_annotations(self) -> None:
        """每个主要科学参数上方都有注释。"""
        src = _get_cell_source(_load_nb(), "03-params-science")
        major_params = [
            "NORMALIZATION_METHOD", "TARGET_SUM", "N_TOP_GENES", "HVG_FLAVOR",
            "BATCH_AWARE_HVG", "EXCLUDE_MT_FROM_HVG", "EXCLUDE_RIBO_FROM_HVG",
            "EXCLUDE_HB_FROM_HVG", "FORCED_INCLUDE_GENES", "REGRESS_OUT",
            "SCALE", "MAX_SCALE_VALUE",
        ]
        for param in major_params:
            # 每个参数赋值前应有四要素注释（含该参数名）
            # 检查参数赋值行之前存在注释行
            lines = src.split("\n")
            assign_line_idx = None
            for i, line in enumerate(lines):
                stripped = line.strip()
                # 匹配 "PARAM = value" 或 "PARAM  = value"（允许 = 前有多个空格）
                if re.match(rf"^{re.escape(param)}\s*=", stripped):
                    assign_line_idx = i
                    break
            if assign_line_idx is None:
                pytest.fail(f"未找到 {param} 的赋值行")
            # 检查赋值行之前的 2-5 行中是否包含注释
            pre_lines = "\n".join(lines[max(0, assign_line_idx - 6):assign_line_idx])
            assert "#" in pre_lines, f"{param} 赋值前缺少注释"
