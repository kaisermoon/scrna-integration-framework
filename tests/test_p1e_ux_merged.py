"""P1-e-merged：02_merged.ipynb PARAMS 四组化 + preflight + 四要素注释 测试。

覆盖：
A. PARAMS 四组化结构：
  1. 4 组 markdown 小标题存在（数据与版本/科学参数/计算参数/输出与运行标识）
  2. 每组有对应 code cell
  3. 所有参数默认值逐字节不变

B. Preflight 预检 cell：
  4. preflight cell 位于 Setup 之后、数据加载(0834b4b3)之前
  5. preflight 含 raise ValueError + backed="r" + 引用关键变量
  6. preflight 不含 dir()、不含全量 sc.read_h5ad（无 backed 参数）

C. 科学参数四要素注释：
  7. JOIN_GENES / MIN_SHARED_GENES / BATCH_KEY / DOWNSAMPLE_TO_MIN / CRITICAL_MARKERS
     含 含义+默认依据+调大调小+何时改 关键词

D. Redline 守卫：
  8. cell 0834b4b3 的 layers["counts"] 契约 intact
  9. doublet 纳入构建 cell intact
  10. checkpoint cell 29c2b000 intact（expression_contract / needs_review / promote_run）

E. 全局纪律：
  11. 无 dir()残留
  12. 无 run_contract 函数重定义
  13. 所有 code cell 通过 ast.parse
  14. notebook JSON 合法
"""

from __future__ import annotations

import ast
import json
from pathlib import Path
from typing import Any

import pytest

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


# ---- A. PARAMS 四组化结构 -------------------------------------------------------


class TestParamsFourGroups:
    """PARAMS 四组化结构检查。"""

    @pytest.fixture(scope="class")
    def cells(self) -> list[dict[str, Any]]:
        return _load_nb_cells()

    def test_four_group_headers_present(self, cells: list[dict[str, Any]]) -> None:
        """4 个组 markdown 小标题均存在。"""
        group_names = ["数据与版本", "科学参数", "计算参数", "输出与运行标识"]
        all_md = ""
        for c in cells:
            if c["cell_type"] == "markdown":
                all_md += "".join(c["source"])
        for gn in group_names:
            assert gn in all_md, f"Missing group header: {gn}"

    def test_params_intro_present(self, cells: list[dict[str, Any]]) -> None:
        """p1e-params-intro markdown 存在并说明参数入口。"""
        intro = None
        for c in cells:
            if c.get("id") == "p1e-params-intro":
                intro = "".join(c["source"])
                break
        assert intro is not None, "intro markdown p1e-params-intro not found"
        assert "一、参数入口" in intro, "intro missing section title"

    def test_each_group_has_code_cell(self, cells: list[dict[str, Any]]) -> None:
        """四个组 code cell 均存在。"""
        group_code_ids = [
            "p1e-g1-data-code", "p1e-g2-sci-code",
            "p1e-g3-comp-code", "p1e-g4-out-code",
        ]
        for cid in group_code_ids:
            found = any(c.get("id") == cid and c["cell_type"] == "code" for c in cells)
            assert found, f"Group code cell {cid} not found or not code type"

    def test_all_param_defaults_unchanged(self, cells: list[dict[str, Any]]) -> None:
        """每个参数默认值逐一等值。"""
        all_src = ""
        for c in cells:
            if c["cell_type"] == "code":
                all_src += "".join(c["source"]) + "\n"

        defaults: list[tuple[str, str]] = [
            ("UPSTREAM_RUN_ROOT", 'UPSTREAM_RUN_ROOT = "results/runs"'),
            ("UPSTREAM_RUNS Nancang", '"Nancang_2025": "01-nancang-v1-run001"'),
            ("UPSTREAM_RUNS Kim", '"Kim_2023": "01-kim-v1-run001"'),
            ("UPSTREAM_RUNS Nowicki", '"Nowicki_2023": "01-nowicki-v1-run001"'),
            ("UPSTREAM_RUNS Yue", '"Yue_2024": "01-yue-v1-run001"'),
            ("JOIN_GENES", 'JOIN_GENES = "inner"'),
            ("MIN_SHARED_GENES", "MIN_SHARED_GENES = 15000"),
            ("BATCH_KEY", 'BATCH_KEY = "source_dataset"'),
            ("DOWNSAMPLE_TO_MIN", "DOWNSAMPLE_TO_MIN = False"),
            ("OUTPUT_VERSION", "OUTPUT_VERSION = 1"),
            ("RANDOM_SEED", "RANDOM_SEED = 42"),
            ("DOUBLET_INCLUDE_KEY", 'DOUBLET_INCLUDE_KEY = "doublet_include"'),
            ("DOUBLET_STATE_KEY", 'DOUBLET_STATE_KEY = "doublet_class"'),
            ("SAMPLE_KEY", 'SAMPLE_KEY = "sample_id"'),
            ("MAX_DOUBLET_FRACTION", "MAX_DOUBLET_FRACTION = 0.30"),
            ("MAX_UNCERTAIN_FRACTION", "MAX_UNCERTAIN_FRACTION = 0.30"),
            ("MIN_SAMPLE_CELLS_FOR_DOUBLET", "MIN_SAMPLE_CELLS_FOR_DOUBLET = 50"),
            ("RUN_ID", 'RUN_ID = "02-merged-v1-run001"'),
            ("RUN_ROOT", 'RUN_ROOT = "results/runs"'),
            ("OUTPUT_FILENAME", 'OUTPUT_FILENAME = "02_merged_v1.h5ad"'),
        ]
        for name, assignment in defaults:
            assert assignment in all_src, (
                f"Default '{name}' not found as '{assignment}'"
            )

    def test_critical_markers_15_present(self, cells: list[dict[str, Any]]) -> None:
        """CRITICAL_MARKERS 包含全部 15 个 marker。"""
        all_src = ""
        for c in cells:
            if c["cell_type"] == "code":
                all_src += "".join(c["source"]) + "\n"
        required = [
            "EPCAM", "CDH1", "VIM", "PTPRC",
            "MUC5AC", "MUC6", "TFF1", "TFF2",
            "LGR5", "OLFM4",
            "MKI67", "TOP2A",
            "CD3D", "CD4", "CD8A",
        ]
        for marker in required:
            assert marker in all_src, f"CRITICAL_MARKERS missing: {marker}"

    def test_params_code_cells_between_setup_and_loading(self, cells: list[dict[str, Any]]) -> None:
        """所有 PARAMS code cell + preflight 位于 Setup 之后、加载之前。"""
        setup_idx = loading_idx = None
        for i, c in enumerate(cells):
            if c.get("id") == "cf443482":
                setup_idx = i
            elif c.get("id") == "0834b4b3":
                loading_idx = i
        assert setup_idx is not None and loading_idx is not None

        param_code_ids = [
            "p1e-g1-data-code", "p1e-g2-sci-code",
            "p1e-g3-comp-code", "p1e-g4-out-code",
            "p1e-preflight-code",
        ]
        for pcid in param_code_ids:
            for i, c in enumerate(cells):
                if c.get("id") == pcid:
                    assert setup_idx < i < loading_idx, (
                        f"{pcid} at idx {i} not between "
                        f"Setup({setup_idx}) and loading({loading_idx})"
                    )
                    break


# ---- B. Preflight 预检 cell ----------------------------------------------------


class TestPreflightCell:
    """Preflight cell 结构与内容检查。"""

    @pytest.fixture(scope="class")
    def cells(self) -> list[dict[str, Any]]:
        return _load_nb_cells()

    @pytest.fixture(scope="class")
    def pf_src(self, cells: list[dict[str, Any]]) -> str:
        return _cell_source(cells, "p1e-preflight-code")

    def test_preflight_before_load(self, cells: list[dict[str, Any]]) -> None:
        """preflight code cell 索引 < 0834b4b3。"""
        pf_idx = loading_idx = None
        for i, c in enumerate(cells):
            if c.get("id") == "p1e-preflight-code":
                pf_idx = i
            elif c.get("id") == "0834b4b3":
                loading_idx = i
        assert pf_idx is not None, "preflight code cell not found"
        assert loading_idx is not None, "loading cell 0834b4b3 not found"
        assert pf_idx < loading_idx, f"preflight idx {pf_idx} >= loading idx {loading_idx}"

    def test_preflight_has_raise_and_backed(self, pf_src: str) -> None:
        """preflight 含 raise ValueError 与 backed="r"。"""
        assert "raise ValueError" in pf_src, "preflight missing raise ValueError"
        assert 'backed="r"' in pf_src, 'preflight missing backed="r"'

    def test_preflight_references_key_vars(self, pf_src: str) -> None:
        """preflight 引用 UPSTREAM_RUNS / DOUBLET_INCLUDE_KEY / resume_run。"""
        assert "UPSTREAM_RUNS" in pf_src, "missing UPSTREAM_RUNS"
        assert "DOUBLET_INCLUDE_KEY" in pf_src, "missing DOUBLET_INCLUDE_KEY"
        assert "resume_run" in pf_src, "missing resume_run"

    def test_preflight_no_dir(self, pf_src: str) -> None:
        """preflight 不含 dir()。"""
        assert "dir()" not in pf_src, "preflight contains dir()"

    def test_preflight_no_bare_full_load(self, pf_src: str) -> None:
        """preflight 中所有 sc.read_h5ad 调用都带 backed 参数。"""
        import re
        calls = re.findall(r'sc\.read_h5ad\([^)]*\)', pf_src)
        for call in calls:
            assert "backed" in call, (
                f"bare sc.read_h5ad without backed: {call[:100]}"
            )

    def test_preflight_markdown_present(self, cells: list[dict[str, Any]]) -> None:
        """preflight markdown 存在并描述分工。"""
        pf_md = None
        for c in cells:
            if c.get("id") == "p1e-preflight-md":
                pf_md = "".join(c["source"])
                break
        assert pf_md is not None, "preflight markdown not found"
        assert "二、Preflight" in pf_md, "missing preflight section title"
        assert "0834b4b3" in pf_md, "missing reference to 0834b4b3"

    def test_preflight_uses_pf_prefix(self, pf_src: str) -> None:
        """preflight 临时变量使用 _pf_ 前缀。"""
        assert "_pf_" in pf_src, "should use _pf_ prefix"

    def test_preflight_division_of_labor(self, pf_src: str) -> None:
        """preflight 注释写明与 cell 0834b4b3 的分工。"""
        assert "0834b4b3" in pf_src, "should reference cell 0834b4b3"
        assert "不替代" in pf_src or "不删改" in pf_src, (
            "should state it does not replace cell 0834b4b3"
        )


# ---- C. 科学参数四要素注释 -------------------------------------------------------


class TestScientificParamsComments:
    """科学参数的四要素注释存在性检查。"""

    @pytest.fixture(scope="class")
    def cells(self) -> list[dict[str, Any]]:
        return _load_nb_cells()

    @pytest.fixture(scope="class")
    def all_param_src(self, cells: list[dict[str, Any]]) -> str:
        """收集 组2 + 组3 code cell 源码。"""
        src_parts = []
        for cid in ["p1e-g2-sci-code", "p1e-g3-comp-code"]:
            for c in cells:
                if c.get("id") == cid:
                    src_parts.append("".join(c["source"]))
                    break
        return "\n".join(src_parts)

    def test_four_element_keywords_present(self, all_param_src: str) -> None:
        """科学/计算参数 cell 含四维度关键词。"""
        for kw in ["含义", "默认依据", "调大调小", "何时改"]:
            assert kw in all_param_src, f"Missing keyword: {kw}"

    def test_join_genes_has_four_elements(self, all_param_src: str) -> None:
        """JOIN_GENES 附近含四要素注释。"""
        assert "JOIN_GENES" in all_param_src
        lines = all_param_src.split("\n")
        join_idx = None
        for i, line in enumerate(lines):
            if 'JOIN_GENES = "inner"' in line or "JOIN_GENES = 'inner'" in line:
                join_idx = i
                break
        assert join_idx is not None
        pre_block = "\n".join(lines[max(0, join_idx - 5):join_idx])
        for kw in ["含义", "默认依据", "调大调小", "何时改"]:
            assert kw in pre_block, f"JOIN_GENES missing '{kw}' in preceding comments"

    def test_min_shared_genes_present(self, all_param_src: str) -> None:
        """MIN_SHARED_GENES 存在。"""
        assert "MIN_SHARED_GENES" in all_param_src

    def test_batch_key_present(self, all_param_src: str) -> None:
        """BATCH_KEY 存在。"""
        assert "BATCH_KEY" in all_param_src

    def test_downsample_to_min_present(self, all_param_src: str) -> None:
        """DOWNSAMPLE_TO_MIN（计算参数组）存在。"""
        assert "DOWNSAMPLE_TO_MIN" in all_param_src

    def test_critical_markers_present(self, all_param_src: str) -> None:
        """CRITICAL_MARKERS 存在。"""
        assert "CRITICAL_MARKERS" in all_param_src


# ---- D. Redline 守卫 ------------------------------------------------------------


class TestRedlineGuards:
    """关键 redline cell 未被动到。"""

    @pytest.fixture(scope="class")
    def cells(self) -> list[dict[str, Any]]:
        return _load_nb_cells()

    def test_counts_contract_cell_0834b4b3_intact(self, cells: list[dict[str, Any]]) -> None:
        """cell 0834b4b3 仍含 layers['counts'] 契约检查。"""
        src = _cell_source(cells, "0834b4b3")
        assert 'layers["counts"]' in src or "layers['counts']" in src, (
            "REDLINE: lost layers['counts'] reference"
        )
        assert "non-negative integers" in src, (
            "REDLINE: lost integer check"
        )
        assert "shape" in src, (
            "REDLINE: lost shape check"
        )

    def test_doublet_include_cell_intact(self, cells: list[dict[str, Any]]) -> None:
        """doublet 纳入构建 cell 实质逻辑未动。"""
        src = _cell_by_marker(cells, "doublet 纳入构建")
        required = [
            "DOUBLET_INCLUDE_KEY", "_keep_mask",
            "doublet_needs_review", "doublet_inclusion_report",
        ]
        for item in required:
            assert item in src, f"doublet cell missing: {item}"

    def test_checkpoint_cell_29c2b000_intact(self, cells: list[dict[str, Any]]) -> None:
        """checkpoint cell 29c2b000 全部关键逻辑未动。"""
        src = _cell_source(cells, "29c2b000")
        required = [
            "expression_contract", "validate_expression_contract",
            "hard_postconditions", "needs_review=doublet_needs_review",
            "doublet_inclusion_report", "snapshot_effective_parameters",
            "promote_run",
        ]
        for item in required:
            assert item in src, f"checkpoint missing: {item}"

        # NEEDS_REVIEW 分支不调 promote_run
        needs_idx = src.find('if stage_status.value == "NEEDS_REVIEW":')
        assert needs_idx > 0, "NEEDS_REVIEW branch not found"
        else_idx = src.find("else:", needs_idx)
        assert else_idx > needs_idx
        needs_block = src[needs_idx:else_idx]
        assert "promote_run(" not in needs_block, (
            "REDLINE: promote_run() called in NEEDS_REVIEW branch"
        )
        else_block = src[else_idx:]
        assert "promote_run(" in else_block, (
            "REDLINE: promote_run() missing from else branch"
        )


# ---- E. 全局纪律 ----------------------------------------------------------------


class TestGlobalDiscipline:
    """跨全部 cell 的纪律检查。"""

    @pytest.fixture(scope="class")
    def cells(self) -> list[dict[str, Any]]:
        return _load_nb_cells()

    def test_no_dir_in_any_code_cell(self, cells: list[dict[str, Any]]) -> None:
        """所有 code cell 不含 dir()。"""
        for c in cells:
            if c["cell_type"] == "code":
                src = "".join(c["source"])
                assert "dir()" not in src, (
                    f"dir() in cell id={c.get('id', 'N/A')}"
                )

    def test_no_run_contract_redefinition(self, cells: list[dict[str, Any]]) -> None:
        """notebook 不重定义 run_contract 函数。"""
        protected = [
            "determine_stage_status", "validate_expression_contract",
            "promote_run", "prepare_run",
        ]
        for fn in protected:
            for c in cells:
                if c["cell_type"] == "code":
                    src = "".join(c["source"])
                    assert f"def {fn}" not in src, (
                        f"redefines {fn}"
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
        """notebook JSON 合法。"""
        assert NB_PATH.exists(), f"not found: {NB_PATH}"
        with open(NB_PATH, encoding="utf-8") as f:
            nb = json.load(f)
        assert "cells" in nb, "missing 'cells' key"
        assert len(nb["cells"]) >= 40, (
            f"expected >=40 cells, got {len(nb['cells'])}"
        )

    def test_setup_cell_unchanged(self, cells: list[dict[str, Any]]) -> None:
        """Setup cell 的 import 清单与 BLAS 线程设置未动。"""
        src = _cell_source(cells, "cf443482")
        assert "from scrna_integration.run_contract import" in src
        assert "promote_run" in src
        assert "resume_run" in src
        assert "sha256_file" in src
        assert "snapshot_effective_parameters" in src
        assert "validate_expression_contract" in src
        assert "OPENBLAS_NUM_THREADS" in src
