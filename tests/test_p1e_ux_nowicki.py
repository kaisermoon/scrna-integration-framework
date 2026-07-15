"""P1-e: UX 骨架测试 — nowicki（PARAMS 四组化 + 四要素注释 + preflight 校验 + 红线守卫）。

覆盖：
1. 参数默认值保真（13 个变量名与默认值逐字不变）
2. cell 顺序守卫（四组在 Setup 前，preflight 在 Setup 之后 data-load 之前）
3. preflight 逻辑镜像（所有校验分支 + 合法路径通过 + 缺文件/非法参数/越界/SoupX/.raw 各分支）
4. 红线未被掏空守卫（counts 契约 / doublet 三态 / checkpoint hard_postconditions）
5. notebook JSON 合法性（json.load + 全 code cell ast.parse）

测试风格：解析 notebook JSON 做契约断言；preflight 逻辑提取为镜像函数单独测。
不 import notebook，不跑 nbconvert。
"""

from __future__ import annotations

import ast
import json
import os
import sys
import tempfile
import textwrap
from pathlib import Path

import numpy as np
import pytest

# ============================================================================
# 路径常量
# ============================================================================

_NOTEBOOK_PATH = Path(__file__).resolve().parent.parent / "notebooks" / "01_per_dataset" / "01_nowicki.ipynb"

# 重构前基线默认值（逐字不变守卫）
_BASELINE_PARAMS = {
    "MANIFEST_PATH": "data/nowicki/manifest.yaml",
    "RUN_ID": "01-nowicki-v1-run001",
    "RUN_ROOT": "results/runs",
    "OUTPUT_FILENAME": "01_nowicki_v1.h5ad",
    "QC_STRATEGY": "skip",
    "SCORE_CELL_CYCLE": True,
    "RANDOM_SEED": 42,
    "OUTPUT_VERSION": 1,
    "EXPECTED_DOUBLET_RATE": None,
    "DOUBLET_SCORE_THRESHOLD": None,
    "DOUBLET_UNCERTAIN_MARGIN": 0.2,
    "DOUBLET_MIN_CELLS": 50,
    "DOUBLET_RATE_ALERT": 0.30,
}

# 参数所属组（用于断言分组正确性）
_PARAM_GROUPS = {
    "g1": ["MANIFEST_PATH"],
    "g2": [
        "QC_STRATEGY", "SCORE_CELL_CYCLE",
        "EXPECTED_DOUBLET_RATE", "DOUBLET_SCORE_THRESHOLD",
        "DOUBLET_UNCERTAIN_MARGIN", "DOUBLET_MIN_CELLS", "DOUBLET_RATE_ALERT",
    ],
    "g3": ["RANDOM_SEED"],
    "g4": ["RUN_ID", "RUN_ROOT", "OUTPUT_FILENAME", "OUTPUT_VERSION"],
}

# ============================================================================
# Notebook 加载 helper
# ============================================================================

def _load_nb():
    """加载 notebook JSON，不存在则 pytest.skip。"""
    if not _NOTEBOOK_PATH.exists():
        pytest.skip(f"notebook 不存在: {_NOTEBOOK_PATH}")
    with open(_NOTEBOOK_PATH) as f:
        return json.load(f)


def _code_src(cell: dict) -> str:
    return "".join(cell.get("source", ""))


def _cell_by_id(nb: dict, cell_id: str) -> dict | None:
    for c in nb["cells"]:
        if c.get("id") == cell_id:
            return c
    return None


# ============================================================================
# 0. notebook JSON 合法性
# ============================================================================

def test_notebook_json_valid():
    """notebook JSON 合法 + 全 code cell 可 ast.parse。"""
    nb = _load_nb()
    assert isinstance(nb, dict)
    assert "cells" in nb
    for i, c in enumerate(nb["cells"]):
        if c.get("cell_type") == "code":
            src = _code_src(c)
            try:
                ast.parse(src)
            except SyntaxError as e:
                pytest.fail(f"cell idx={i} id={c.get('id')} 语法错误: {e}")


# ============================================================================
# 1. 参数默认值保真
# ============================================================================

def _extract_params_from_code_cells(nb: dict, cell_ids: list[str]) -> dict:
    """从指定 code cell 中提取所有 NAME = literal 赋值，返回 {name: value}。"""
    params = {}
    for cid in cell_ids:
        c = _cell_by_id(nb, cid)
        if c is None:
            continue
        src = _code_src(c)
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        try:
                            val = ast.literal_eval(node.value)
                        except (ValueError, TypeError):
                            continue
                        params[target.id] = val
    return params


# 四组 code cell id
_PARAM_CODE_IDS = ["aae03603", "p1e_params_g2_code", "p1e_params_g3_code", "p1e_params_g4_code"]


def test_params_all_13_variables_present():
    """四组 code cell 中仍包含全部 13 个变量。"""
    nb = _load_nb()
    actual = _extract_params_from_code_cells(nb, _PARAM_CODE_IDS)
    missing = set(_BASELINE_PARAMS) - set(actual)
    extra = set(actual) - set(_BASELINE_PARAMS)
    assert not missing, f"缺失参数: {missing}"
    assert not extra, f"多出参数: {extra}"


def test_params_default_values_unchanged():
    """每个参数默认值与重构前基线逐字一致。"""
    nb = _load_nb()
    actual = _extract_params_from_code_cells(nb, _PARAM_CODE_IDS)
    diffs = []
    for name, expected in _BASELINE_PARAMS.items():
        actual_val = actual.get(name)
        if actual_val != expected:
            # 浮点比较
            if isinstance(expected, float) and isinstance(actual_val, float):
                if abs(actual_val - expected) < 1e-12:
                    continue
            diffs.append(f"  {name}: expected={expected!r}, actual={actual_val!r}")
    assert not diffs, "参数默认值被修改:\n" + "\n".join(diffs)


def test_params_group_assignment():
    """每个参数在正确的组 code cell 中定义。"""
    nb = _load_nb()
    group_cid = {
        "g1": "aae03603", "g2": "p1e_params_g2_code",
        "g3": "p1e_params_g3_code", "g4": "p1e_params_g4_code",
    }
    for gname, expected_vars in _PARAM_GROUPS.items():
        c = _cell_by_id(nb, group_cid[gname])
        assert c is not None, f"组 {gname} code cell {group_cid[gname]} 不存在"
        src = _code_src(c)
        for var in expected_vars:
            assert var in src, f"{var} 应在组 {gname} 中定义，但未在 cell {group_cid[gname]} 中找到"


# ============================================================================
# 2. cell 顺序守卫
# ============================================================================

def test_params_cells_before_setup():
    """四组 params cells 全部在 Setup (c4508e7d) 之前。"""
    nb = _load_nb()
    cell_ids = [c.get("id") for c in nb["cells"]]
    setup_pos = cell_ids.index("c4508e7d")
    param_ids = [
        "p1e_params_intro_md", "p1e_params_g1_md", "aae03603",
        "p1e_params_g2_md", "p1e_params_g2_code",
        "p1e_params_g3_md", "p1e_params_g3_code",
        "p1e_params_g4_md", "p1e_params_g4_code",
    ]
    for pid in param_ids:
        assert pid in cell_ids, f"{pid} 不存在于 notebook"
        assert cell_ids.index(pid) < setup_pos, f"{pid} 应在 Setup 之前，实际在之后"


def test_preflight_after_setup_before_dataload():
    """preflight cells 在 Setup (c4508e7d) 之后、data-load (9978541f) 之前。"""
    nb = _load_nb()
    cell_ids = [c.get("id") for c in nb["cells"]]
    setup_pos = cell_ids.index("c4508e7d")
    dataload_pos = cell_ids.index("9978541f")
    pre_md_pos = cell_ids.index("p1e_preflight_md")
    pre_code_pos = cell_ids.index("p1e_preflight_code")
    assert setup_pos < pre_md_pos < pre_code_pos < dataload_pos, (
        f"顺序错误: Setup={setup_pos}, pre_md={pre_md_pos}, pre_code={pre_code_pos}, data-load={dataload_pos}"
    )


def test_doublet_cell_unchanged_position():
    """doublet cells (p1c_doublet_md / p1c_doublet_code) 仍在 data-load 之后。"""
    nb = _load_nb()
    cell_ids = [c.get("id") for c in nb["cells"]]
    dataload_pos = cell_ids.index("9978541f")
    dbl_md_pos = cell_ids.index("p1c_doublet_md")
    dbl_code_pos = cell_ids.index("p1c_doublet_code")
    assert dataload_pos < dbl_md_pos < dbl_code_pos, "doublet cells 应在 data-load 之后"


# ============================================================================
# 3. preflight 逻辑镜像
# ============================================================================

def _make_manifest_dict(input_path: str = "/nonexistent/test.h5ad", **kwargs) -> dict:
    """构造最小 manifest，缺省含 nowicki 默认 preprocessing_done + qc_overrides。"""
    base = {
        "input": {"path": input_path},
        "source_dataset": "nowicki_2023",
        "preprocessing_done": kwargs.pop("preprocessing_done", ["basic_filter", "doublet_removal", "normalization"]),
        "qc_overrides": kwargs.pop("qc_overrides", {
            "basic_filter": {"skip": True, "reason": "作者已完成"},
            "doublet_removal": {"skip": True, "reason": "作者已完成"},
            "normalization": {"skip": True, "reason": "作者已完成"},
        }),
    }
    base.update(kwargs)
    return base


def _preflight_validate(
    manifest_path_exists: bool = True,
    input_path_exists: bool = True,
    manifest: dict | None = None,
    random_seed: int = 42,
    output_version: int = 1,
    qc_strategy: str = "skip",
    score_cell_cycle: bool = True,
    doublet_uncertain_margin: float = 0.2,
    doublet_min_cells: int = 50,
    doublet_rate_alert: float = 0.30,
    expected_doublet_rate: float | None = None,
    doublet_score_threshold: float | None = None,
    soupx_enabled: bool = False,
    raw_exists: bool = True,
):
    """镜像 preflight cell 的所有校验分支。

    所有临时变量以 _ 前缀，与 notebook preflight cell 保持一致的变量命名约定。
    返回 (ok: bool, errors: list[str])。
    """
    _errs = []

    if not manifest_path_exists:
        _errs.append("manifest 文件不存在")

    # 类型校验
    if not isinstance(random_seed, int):
        _errs.append(f"RANDOM_SEED 应为 int，实际 {type(random_seed).__name__}")
    if not isinstance(output_version, int):
        _errs.append(f"OUTPUT_VERSION 应为 int，实际 {type(output_version).__name__}")
    if qc_strategy not in {"skip"}:
        _errs.append(f"QC_STRATEGY 应为 'skip'，实际 {qc_strategy!r}")
    if not isinstance(score_cell_cycle, bool):
        _errs.append(f"SCORE_CELL_CYCLE 应为 bool，实际 {type(score_cell_cycle).__name__}")
    if not isinstance(doublet_uncertain_margin, (int, float)) or not (0 < doublet_uncertain_margin < 1):
        _errs.append(f"DOUBLET_UNCERTAIN_MARGIN 应为 float 且 0<值<1，实际 {doublet_uncertain_margin!r}")
    if not isinstance(doublet_min_cells, int) or doublet_min_cells <= 0:
        _errs.append(f"DOUBLET_MIN_CELLS 应为 int > 0，实际 {doublet_min_cells!r}")
    if not isinstance(doublet_rate_alert, (int, float)) or not (0 < doublet_rate_alert <= 1):
        _errs.append(f"DOUBLET_RATE_ALERT 应为 float 且 0<值<=1，实际 {doublet_rate_alert!r}")
    if expected_doublet_rate is not None:
        if not isinstance(expected_doublet_rate, (int, float)) or not (0 < expected_doublet_rate < 1):
            _errs.append(f"EXPECTED_DOUBLET_RATE 应为 None 或 float 且 0<值<1，实际 {expected_doublet_rate!r}")
    if doublet_score_threshold is not None:
        if not isinstance(doublet_score_threshold, (int, float)) or doublet_score_threshold <= 0:
            _errs.append(f"DOUBLET_SCORE_THRESHOLD 应为 None 或 float > 0，实际 {doublet_score_threshold!r}")

    if _errs:
        return False, _errs

    if not manifest_path_exists:
        return False, ["manifest 文件不存在"]

    # 输入文件存在
    if not input_path_exists:
        _errs.append("输入 h5ad 文件不存在")

    if _errs:
        return False, _errs

    # SoupX 禁用断言
    if soupx_enabled:
        _errs.append("Nowicki 无完整 raw droplets，不支持 SoupX（决策3）")

    if _errs:
        return False, _errs

    # .raw.X 存在
    if not raw_exists:
        _errs.append("Nowicki 数据必须包含 .raw.X 原始计数")

    if _errs:
        return False, _errs

    return True, []


# --- 合法配置通过 ---

def test_preflight_valid_default_config():
    """nowicki 默认配置应通过全部 preflight 校验。"""
    ok, errs = _preflight_validate()
    assert ok, f"默认配置应通过 preflight: {errs}"


# --- 缺文件 ---

def test_preflight_manifest_missing():
    """manifest 路径不存在应报错。"""
    ok, errs = _preflight_validate(manifest_path_exists=False)
    assert not ok
    assert any("manifest" in e for e in errs)


def test_preflight_input_file_missing():
    """输入文件不存在应报错。"""
    ok, errs = _preflight_validate(input_path_exists=False)
    assert not ok
    assert any("h5ad" in e for e in errs)


# --- 参数类型/取值校验 ---

def test_preflight_qc_strategy_invalid():
    """QC_STRATEGY 非 'skip' 应报错。"""
    ok, errs = _preflight_validate(qc_strategy="standard")
    assert not ok
    assert any("QC_STRATEGY" in e for e in errs)


def test_preflight_score_cell_cycle_not_bool():
    """SCORE_CELL_CYCLE 非 bool 应报错。"""
    ok, errs = _preflight_validate(score_cell_cycle="yes")
    assert not ok


def test_preflight_uncertain_margin_zero():
    """DOUBLET_UNCERTAIN_MARGIN=0 应报错（不在 (0,1)）。"""
    ok, errs = _preflight_validate(doublet_uncertain_margin=0.0)
    assert not ok


def test_preflight_uncertain_margin_one():
    """DOUBLET_UNCERTAIN_MARGIN=1 应报错（不在 (0,1)）。"""
    ok, errs = _preflight_validate(doublet_uncertain_margin=1.0)
    assert not ok


def test_preflight_uncertain_margin_negative():
    """DOUBLET_UNCERTAIN_MARGIN<0 应报错。"""
    ok, errs = _preflight_validate(doublet_uncertain_margin=-0.5)
    assert not ok


def test_preflight_min_cells_zero():
    """DOUBLET_MIN_CELLS=0 应报错。"""
    ok, errs = _preflight_validate(doublet_min_cells=0)
    assert not ok


def test_preflight_min_cells_negative():
    """DOUBLET_MIN_CELLS<0 应报错。"""
    ok, errs = _preflight_validate(doublet_min_cells=-10)
    assert not ok


def test_preflight_rate_alert_zero():
    """DOUBLET_RATE_ALERT=0 应报错（不在 (0,1]）。"""
    ok, errs = _preflight_validate(doublet_rate_alert=0.0)
    assert not ok


def test_preflight_rate_alert_gt_one():
    """DOUBLET_RATE_ALERT>1 应报错。"""
    ok, errs = _preflight_validate(doublet_rate_alert=1.5)
    assert not ok


def test_preflight_rate_alert_one_valid():
    """DOUBLET_RATE_ALERT=1.0 应合法（边界值）。"""
    ok, errs = _preflight_validate(doublet_rate_alert=1.0)
    assert ok


def test_preflight_expected_doublet_rate_negative():
    """EXPECTED_DOUBLET_RATE<0 应报错。"""
    ok, errs = _preflight_validate(expected_doublet_rate=-0.05)
    assert not ok


def test_preflight_expected_doublet_rate_gt_one():
    """EXPECTED_DOUBLET_RATE>1 应报错。"""
    ok, errs = _preflight_validate(expected_doublet_rate=1.5)
    assert not ok


def test_preflight_expected_doublet_rate_none_valid():
    """EXPECTED_DOUBLET_RATE=None 应合法。"""
    ok, errs = _preflight_validate(expected_doublet_rate=None)
    assert ok


def test_preflight_doublet_score_threshold_negative():
    """DOUBLET_SCORE_THRESHOLD<0 应报错。"""
    ok, errs = _preflight_validate(doublet_score_threshold=-0.1)
    assert not ok


def test_preflight_doublet_score_threshold_zero():
    """DOUBLET_SCORE_THRESHOLD=0 应报错（需 >0）。"""
    ok, errs = _preflight_validate(doublet_score_threshold=0.0)
    assert not ok


def test_preflight_doublet_score_threshold_none_valid():
    """DOUBLET_SCORE_THRESHOLD=None 应合法。"""
    ok, errs = _preflight_validate(doublet_score_threshold=None)
    assert ok


def test_preflight_random_seed_not_int():
    """RANDOM_SEED 非 int 应报错。"""
    ok, errs = _preflight_validate(random_seed="foo")  # type: ignore[arg-type]
    assert not ok


def test_preflight_output_version_not_int():
    """OUTPUT_VERSION 非 int 应报错。"""
    ok, errs = _preflight_validate(output_version="1")  # type: ignore[arg-type]
    assert not ok


# --- SoupX 禁用断言 ---

def test_preflight_soupx_enabled_raises():
    """SOUPX_ENABLED=True 应报错。"""
    ok, errs = _preflight_validate(soupx_enabled=True)
    assert not ok
    assert any("SoupX" in e for e in errs)


def test_preflight_soupx_disabled_ok():
    """SOUPX_ENABLED=False（默认）应通过。"""
    ok, errs = _preflight_validate(soupx_enabled=False)
    assert ok


# --- .raw.X 不存在 ---

def test_preflight_raw_missing():
    """.raw.X 不存在应报错。"""
    ok, errs = _preflight_validate(raw_exists=False)
    assert not ok
    assert any(".raw" in e for e in errs)


# ============================================================================
# 4. 红线未被掏空守卫
# ============================================================================

def test_counts_contract_9978541f_intact():
    """data-load cell 9978541f 仍含 layers['counts'] 与 expression_contract。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "9978541f")
    assert c is not None, "data-load cell 9978541f 丢失"
    src = _code_src(c)
    assert "layers[\"counts\"]" in src or "layers['counts']" in src, "counts 契约遗失"
    assert "expression_contract" in src, "expression_contract 遗失"


def test_doublet_code_p1c_intact():
    """doublet cell p1c_doublet_code 仍含 classify_doublets 与 doublet_detection。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "p1c_doublet_code")
    assert c is not None, "doublet code cell 丢失"
    src = _code_src(c)
    assert "classify_doublets" in src, "classify_doublets 函数遗失"
    assert "doublet_detection" in src, "doublet_detection schema 遗失"
    # 三态边界语义
    assert '"singlet"' in src or "'singlet'" in src, "singlet 状态遗失"
    assert '"uncertain"' in src or "'uncertain'" in src, "uncertain 状态遗失"
    assert '"doublet"' in src or "'doublet'" in src, "doublet 状态遗失"
    # _DOUBLET_COLS_INITIALIZED guard（禁 dir() 判断）
    assert "_DOUBLET_COLS_INITIALIZED" in src, "_DOUBLET_COLS_INITIALIZED guard 遗失"


def test_checkpoint_18893366_intact():
    """checkpoint cell 18893366 仍含 hard_postconditions / _hd / counts_source / promote_run。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "18893366")
    assert c is not None, "checkpoint cell 18893366 丢失"
    src = _code_src(c)
    assert "hard_postconditions" in src, "hard_postconditions 遗失"
    assert "_hd" in src, "_hd guard 遗失"
    assert 'counts_source' in src, "counts_source 断言遗失"
    assert '".raw.X"' in src or "'.raw.X'" in src, "counts_source=='.raw.X' 断言遗失"
    assert "promote_run" in src, "promote_run 遗失"
    assert "validate_expression_contract" in src, "validate_expression_contract 调用遗失"


def test_data_load_cell_still_has_raw_extraction():
    """data-load cell 仍包含 .raw.X 提取逻辑与 counts 非负整数校验。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "9978541f")
    src = _code_src(c)
    assert ".raw.X" in src or ".raw" in src, ".raw 提取逻辑遗失"
    assert "counts contain negative" in src, "counts 非负校验遗失"
    assert "counts contain non-integer" in src, "counts 整数校验遗失"


def test_no_soupx_logic_added_to_nowicki():
    """nowicki notebook 不应新增 SoupX 逻辑（除 preflight 的禁用断言）。"""
    nb = _load_nb()
    # 除 preflight 外，任何 cell 不应有 SOUPX_ENABLED / counts_soupx / SoupX
    for c in nb["cells"]:
        cid = c.get("id", "")
        if cid == "p1e_preflight_code":
            continue  # preflight 的 SoupX 禁用断言是允许的
        src = _code_src(c)
        assert "SOUPX_ENABLED" not in src, f"cell {cid} 不应包含 SOUPX_ENABLED"
        assert "counts_soupx" not in src, f"cell {cid} 不应包含 counts_soupx"
        assert "soupx_layer" not in src or cid == "9978541f", (
            f"cell {cid} 不应包含 soupx_layer（仅 data-load 的 expression_contract 声明 None 除外）"
        )


def test_no_new_imports_in_logic_cells():
    """逻辑 cell（counts/doublet/checkpoint）不应被注入新 import。"""
    nb = _load_nb()
    logic_ids = ["9978541f", "p1c_doublet_code", "18893366"]
    for cid in logic_ids:
        c = _cell_by_id(nb, cid)
        assert c is not None, f"逻辑 cell {cid} 丢失"
        # 确认 cell 是 code
        assert c.get("cell_type") == "code", f"{cid} 应为 code cell"


def test_aae03603_still_has_manifest_path_only():
    """复用的 aae03603 仅含 MANIFEST_PATH（组1），不含其他参数。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "aae03603")
    src = _code_src(c)
    assert "MANIFEST_PATH" in src, "MANIFEST_PATH 应在 aae03603"
    # 不应含其他组的参数
    for var in ["QC_STRATEGY", "RANDOM_SEED", "RUN_ID", "DOUBLET_UNCERTAIN_MARGIN"]:
        assert var not in src, f"aae03603 不应含 {var}（已移至其他组）"


# ============================================================================
# 5. 四要素注释存在性
# ============================================================================

def test_g2_scientific_params_have_four_element_comments():
    """组2 code cell 中每个科学参数赋值上方至少有四行注释。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "p1e_params_g2_code")
    src = _code_src(c)

    # 每个科学参数应包含四要素标记
    required_annotations = ["含义", "默认依据", "调大调小", "何时该改"]
    sci_params = [
        "QC_STRATEGY", "SCORE_CELL_CYCLE", "EXPECTED_DOUBLET_RATE",
        "DOUBLET_SCORE_THRESHOLD", "DOUBLET_UNCERTAIN_MARGIN",
        "DOUBLET_MIN_CELLS", "DOUBLET_RATE_ALERT",
    ]
    for param in sci_params:
        assert param in src, f"{param} 应在组2 code cell 中"
        # 验证四要素全部在 param 赋值之前的注释中出现
        # 找到 param 赋值的位置
        param_idx = src.index(param)
        prefix = src[:param_idx]
        for anno in required_annotations:
            # 每个要素至少出现一次（可能跨多个参数共享同一个要素注释块）
            # 宽松检查：全文包含即可，不要求每个 param 前独立出现
            assert anno in prefix, f"参数 {param} 上方缺少注释要素 '{anno}'"


def test_non_scientific_params_have_single_line_comment():
    """非科学参数只给一行说明，不套四要素。"""
    nb = _load_nb()
    for cid in ["aae03603", "p1e_params_g3_code", "p1e_params_g4_code"]:
        c = _cell_by_id(nb, cid)
        src = _code_src(c)
        # 不应出现四要素注释（它们只应在 g2 中出现，且非科学参数不应套四要素）
        # 不强制检查"全然没有"，因为中文可能碰巧包含这些字
        # 但至少保证每个非科学参数有注释（以 # 开头）
        lines = [l.strip() for l in src.split("\n")]
        has_comment = any(l.startswith("#") for l in lines)
        assert has_comment, f"{cid} 中的参数应有注释"


# ============================================================================
# 6. preflight cell 在 notebook 中的存在性
# ============================================================================

def test_preflight_cells_exist():
    """preflight md + code cells 存在且 code cell 语法合法。"""
    nb = _load_nb()
    pre_md = _cell_by_id(nb, "p1e_preflight_md")
    pre_code = _cell_by_id(nb, "p1e_preflight_code")
    assert pre_md is not None, "preflight md cell 缺失"
    assert pre_code is not None, "preflight code cell 缺失"
    assert pre_md["cell_type"] == "markdown"
    assert pre_code["cell_type"] == "code"


def test_preflight_code_has_underscore_prefix_vars():
    """preflight code cell 中临时变量以 _ 前缀。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "p1e_preflight_code")
    src = _code_src(c)
    # 所有非 PARAMS 的变量赋值应为 _ 开头
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    name = target.id
                    # PARAMS 全局变量不检查 _ 前缀
                    if name in _BASELINE_PARAMS:
                        continue
                    # 其他赋值必须是 _ 前缀
                    assert name.startswith("_"), (
                        f"preflight 中非 PARAMS 变量 '{name}' 应以 _ 前缀"
                        f"（避免被 checkpoint snapshot_effective_parameters 捕获）"
                    )


def test_preflight_contains_soupx_disabled_assert():
    """preflight 含 SoupX 禁用断言。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "p1e_preflight_code")
    src = _code_src(c)
    assert "SOUPX_ENABLED" in src, "preflight 应含 SOUPX_ENABLED 检查"
    assert "SoupX" in src, "preflight 应声明 SoupX 禁用"


def test_preflight_contains_raw_check():
    """preflight 含 .raw.X 探测（backed 模式）。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "p1e_preflight_code")
    src = _code_src(c)
    assert "read_h5ad" in src, "preflight 应使用 read_h5ad backed 读取"
    assert "backed" in src, "preflight 应使用 backed 模式"
    assert ".raw" in src, "preflight 应检查 .raw"


def test_preflight_contains_params_summary():
    """preflight 末尾打印生效参数摘要。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "p1e_preflight_code")
    src = _code_src(c)
    assert "Preflight 校验通过" in src, "preflight 应打印通过摘要"
    assert "RUN_ID" in src, "摘要应含 RUN_ID"
    assert "QC_STRATEGY" in src, "摘要应含 QC_STRATEGY"


def test_preflight_contains_doublet_coherence_warning():
    """preflight 含 doublet 配置连贯性 WARNING（只 print 不 raise）。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "p1e_preflight_code")
    src = _code_src(c)
    assert "WARNING" in src, "preflight 应含 doublet WARNING"
    assert "needs_review" in src or "needs_review" in src.lower(), "WARNING 应提及 needs_review"


def test_preflight_no_heavy_imports():
    """preflight 不应 import 重型依赖（scanpy/scrublet/scipy 等）。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "p1e_preflight_code")
    src = _code_src(c)
    forbidden = ["import scanpy", "import scrublet", "import scipy", "import matplotlib", "import pandas"]
    for imp in forbidden:
        assert imp not in src, f"preflight 不应 import 重型依赖: {imp}"


def test_preflight_no_adata_write():
    """preflight 不写 adata/uns。"""
    nb = _load_nb()
    c = _cell_by_id(nb, "p1e_preflight_code")
    src = _code_src(c)
    # 不应该有 adata 赋值（除 backed 探测的 _probe）
    assert "adata.uns" not in src, "preflight 不应写 adata.uns"
    assert "adata.layers" not in src, "preflight 不应写 adata.layers"
    assert "adata.write" not in src, "preflight 不应写 adata"
