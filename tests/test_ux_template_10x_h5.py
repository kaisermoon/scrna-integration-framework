"""测试 10x h5 格式模板：01_template_10x_h5.ipynb 的 UX guided 骨架收尾。

覆盖 PARAMS 四组化、preflight cell、科学参数四要素注释与红线守护。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# helpers（参照 test_doublet_10x_h5.py）
# ---------------------------------------------------------------------------

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


def _find_cell_idx(nb: dict, marker: str, cell_type: str = "code") -> int:
    """按内容 marker 定位 cell，返回 index。"""
    for i, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != cell_type:
            continue
        src = "".join(cell["source"])
        if marker in src:
            return i
    raise LookupError(f"cell with marker {marker!r} not found in notebook")


# ---------------------------------------------------------------------------
# A. PARAMS 回归守护
# ---------------------------------------------------------------------------

class TestParamsFourGroupRegression:
    """验证 PARAMS cell 四组化后所有变量仍在，且默认值等于现值。"""

    def test_params_all_variables_present(self):
        """exec PARAMS cell 后四组全部变量都在 env 中。"""
        nb = _load_notebook()
        env = _cell(nb, "=== PARAMS ===")
        # 组一
        assert isinstance(env.get("MANIFEST_PATH"), str) and env["MANIFEST_PATH"].endswith((".yaml", ".yml"))
        # 组二
        assert env.get("QC_STRATEGY") == "adaptive"
        assert env.get("N_MAD") == 3
        assert env.get("PER_SAMPLE_MAD") is True
        assert env.get("MIN_CELLS_PER_GENE") == 3
        assert "DOUBLET_SCORE_THRESHOLD" in env
        assert "DOUBLET_SCORE_HIGH" in env
        assert "DOUBLET_SCORE_LOW" in env
        assert env.get("DOUBLET_RATE_ALERT_HIGH") == 0.40
        assert "DOUBLET_RATE_ALERT_LOW" in env
        assert env.get("DOUBLET_MIN_CELLS") == 50
        # 固定阈值占位
        assert env.get("MIN_GENES") == 200
        assert env.get("MAX_GENES") == 6000
        assert env.get("MIN_COUNTS") == 500
        assert env.get("MAX_PCT_MT") == 20
        # 组三
        assert env.get("SOUPX_ENABLED") is False
        assert "EXPECTED_DOUBLET_RATE" in env
        assert env.get("DOUBLET_UNCERTAIN_INCLUDE") is True
        assert env.get("SCORE_CELL_CYCLE") is True
        assert env.get("FLAG_HEMOGLOBIN") is False
        assert env.get("FLAG_STRESS_GENES") is True
        # 组四
        assert isinstance(env.get("RUN_ID"), str) and env["RUN_ID"].startswith("01-")
        assert env.get("RUN_ROOT") == "results/runs"
        assert isinstance(env.get("OUTPUT_FILENAME"), str) and env["OUTPUT_FILENAME"].startswith("01_") and env["OUTPUT_FILENAME"].endswith(".h5ad")
        assert env.get("OUTPUT_VERSION") == 1
        assert env.get("RANDOM_SEED") == 42

    def test_params_defaults_unchanged(self):
        """逐个变量断言默认值与原 notebook 一致——四组化只重排不改变语义。"""
        nb = _load_notebook()
        env = _cell(nb, "=== PARAMS ===")
        defaults = {
            "QC_STRATEGY": "adaptive",
            "N_MAD": 3,
            "PER_SAMPLE_MAD": True,
            "MIN_CELLS_PER_GENE": 3,
            "SOUPX_ENABLED": False,
            "EXPECTED_DOUBLET_RATE": None,
            "DOUBLET_UNCERTAIN_INCLUDE": True,
            "SCORE_CELL_CYCLE": True,
            "FLAG_HEMOGLOBIN": False,
            "FLAG_STRESS_GENES": True,
            "DOUBLET_SCORE_THRESHOLD": None,
            "DOUBLET_SCORE_HIGH": None,
            "DOUBLET_SCORE_LOW": None,
            "DOUBLET_RATE_ALERT_HIGH": 0.40,
            "DOUBLET_RATE_ALERT_LOW": None,
            "DOUBLET_MIN_CELLS": 50,
            "OUTPUT_VERSION": 1,
            "RANDOM_SEED": 42,
            "MIN_GENES": 200,
            "MAX_GENES": 6000,
            "MIN_COUNTS": 500,
            "MAX_PCT_MT": 20,
        }
        for key, expected in defaults.items():
            actual = env[key]
            assert actual == expected, (
                f"{key}: expected {expected!r}, got {actual!r}"
            )

    def test_no_variables_lost(self):
        """四组化后 PARAMS cell exec 的变量数量不少于关键变量全集。"""
        nb = _load_notebook()
        env = _cell(nb, "=== PARAMS ===")
        required_keys = {
            "MANIFEST_PATH",
            "QC_STRATEGY", "N_MAD", "PER_SAMPLE_MAD", "MIN_CELLS_PER_GENE",
            "DOUBLET_SCORE_THRESHOLD", "DOUBLET_SCORE_HIGH", "DOUBLET_SCORE_LOW",
            "DOUBLET_RATE_ALERT_HIGH", "DOUBLET_RATE_ALERT_LOW", "DOUBLET_MIN_CELLS",
            "MIN_GENES", "MAX_GENES", "MIN_COUNTS", "MAX_PCT_MT",
            "SOUPX_ENABLED", "EXPECTED_DOUBLET_RATE", "DOUBLET_UNCERTAIN_INCLUDE",
            "SCORE_CELL_CYCLE", "FLAG_HEMOGLOBIN", "FLAG_STRESS_GENES",
            "RUN_ID", "RUN_ROOT", "OUTPUT_FILENAME", "OUTPUT_VERSION", "RANDOM_SEED",
        }
        missing = required_keys - set(env.keys())
        assert not missing, f"PARAMS cell missing keys: {missing}"


# ---------------------------------------------------------------------------
# B. 四组 banner 存在
# ---------------------------------------------------------------------------

class TestFourGroupBanners:
    """验证 PARAMS cell 源码中含四条四组 banner 注释。"""

    def test_four_banners_present(self):
        nb = _load_notebook()
        src = _source(nb, "=== PARAMS ===")
        banners = [
            "组一 · 数据源与版本",
            "组二 · QC 与科学阈值",
            "组三 · 方法开关",
            "组四 · 输出与运行标识",
        ]
        for b in banners:
            assert b in src, f"Banner missing: {b}"


# ---------------------------------------------------------------------------
# C. params-guide markdown 存在
# ---------------------------------------------------------------------------

class TestParamsGuideMarkdown:
    """验证「唯一参数入口」说明 markdown cell 存在。"""

    def test_params_guide_cell_exists(self):
        nb = _load_notebook()
        cell = _find_cell_by_marker(nb, "唯一参数入口", cell_type="markdown")
        assert cell["id"] == "p1e-params-guide", f"unexpected id: {cell['id']}"
        src = "".join(cell["source"])
        assert "四组" in src
        assert "数据源与版本" in src
        assert "QC 与科学阈值" in src
        assert "方法开关" in src
        assert "输出与运行标识" in src
        assert "RUN_ID" in src

    def test_params_guide_before_params_cell(self):
        """params-guide markdown 在 PARAMS code cell 之前。"""
        nb = _load_notebook()
        guide_idx = _find_cell_idx(nb, "唯一参数入口", cell_type="markdown")
        params_idx = _find_cell_idx(nb, "=== PARAMS ===", cell_type="code")
        assert guide_idx < params_idx, (
            f"params-guide ({guide_idx}) should be before PARAMS cell ({params_idx})"
        )


# ---------------------------------------------------------------------------
# D. preflight cell 存在且顺序正确
# ---------------------------------------------------------------------------

class TestPreflightCellExistsAndOrder:
    """验证 preflight cell 存在且位于 bootstrap 与 data-load 之间。"""

    def test_preflight_code_cell_exists(self):
        nb = _load_notebook()
        cell = _find_cell_by_marker(nb, "Preflight：执行前校验", cell_type="code")
        assert cell["id"] == "p1e-preflight", f"unexpected id: {cell['id']}"

    def test_preflight_md_cell_exists(self):
        nb = _load_notebook()
        cell = _find_cell_by_marker(nb, "Preflight 执行前校验", cell_type="markdown")
        assert cell["id"] == "p1e-preflight-md", f"unexpected id: {cell['id']}"

    def test_preflight_after_bootstrap(self):
        nb = _load_notebook()
        bootstrap_idx = _find_cell_idx(nb, "启动脚手架")
        preflight_idx = _find_cell_idx(nb, "Preflight：执行前校验")
        assert preflight_idx > bootstrap_idx, (
            f"preflight ({preflight_idx}) should be after bootstrap ({bootstrap_idx})"
        )

    def test_preflight_before_dataload(self):
        nb = _load_notebook()
        preflight_idx = _find_cell_idx(nb, "Preflight：执行前校验")
        dataload_idx = _find_cell_idx(nb, "数据读入")
        assert preflight_idx < dataload_idx, (
            f"preflight ({preflight_idx}) should be before data-load ({dataload_idx})"
        )


# ---------------------------------------------------------------------------
# E. preflight 合法路径通过
# ---------------------------------------------------------------------------

class TestPreflightHappyPath:
    """用 tmp_path 造 manifest + dummy h5ad，exec preflight 不抛错。"""

    def test_preflight_passes_with_valid_inputs(self, tmp_path):
        """构造合法 manifest 与空 h5ad 文件，exec preflight 应打印通过信息。"""
        import yaml

        # 创建 dummy h5ad 文件（空字节即可，preflight 只检查 is_file 不读内容）
        dummy_h5ad = tmp_path / "dummy.h5ad"
        dummy_h5ad.write_bytes(b"")

        # 创建合法 manifest
        manifest_content = {
            "input": {"format": "10x_h5", "path": str(dummy_h5ad)},
            "source_dataset": "kim_test",
        }
        manifest_path = tmp_path / "manifest.yaml"
        with open(manifest_path, "w") as f:
            yaml.dump(manifest_content, f)

        # 构造 PARAMS env（用默认值，仅改 MANIFEST_PATH）
        nb = _load_notebook()
        params_env = _cell(nb, "=== PARAMS ===")
        params_env["MANIFEST_PATH"] = str(manifest_path)
        # EXPECTED_DOUBLET_RATE 默认为 None（不触发检测），符合 preflight 预期
        # SOUPX_ENABLED 默认 False，符合 Kim 数据集特点

        # 提供 Path 与 yaml（模拟 bootstrap 已导入）
        params_env["yaml"] = yaml
        params_env["Path"] = Path

        # exec preflight cell
        src = _source(nb, "Preflight：执行前校验")
        exec(src, params_env)
        # 不抛错即通过


# ---------------------------------------------------------------------------
# F. preflight 非法即 raise（逐分支）
# ---------------------------------------------------------------------------

class TestPreflightErrorPaths:
    """逐分支验证 preflight 在非法输入下 raise 正确异常。"""

    def test_manifest_not_found(self):
        """MANIFEST_PATH 指向不存在文件 → FileNotFoundError。"""
        nb = _load_notebook()
        params_env = _cell(nb, "=== PARAMS ===")
        params_env["MANIFEST_PATH"] = "/nonexistent/path/manifest.yaml"
        params_env["yaml"] = __import__("yaml")
        params_env["Path"] = Path

        src = _source(nb, "Preflight：执行前校验")
        with pytest.raises(FileNotFoundError, match="manifest"):
            exec(src, params_env)

    def test_qc_strategy_invalid(self, tmp_path):
        """QC_STRATEGY='bad' → ValueError。"""
        import yaml

        manifest_path = tmp_path / "manifest.yaml"
        dummy_h5ad = tmp_path / "dummy.h5ad"
        dummy_h5ad.write_bytes(b"")
        with open(manifest_path, "w") as f:
            yaml.dump({"input": {"format": "10x_h5", "path": str(dummy_h5ad)}, "source_dataset": "test"}, f)

        nb = _load_notebook()
        params_env = _cell(nb, "=== PARAMS ===")
        params_env["MANIFEST_PATH"] = str(manifest_path)
        params_env["QC_STRATEGY"] = "bad"
        params_env["yaml"] = yaml
        params_env["Path"] = Path

        src = _source(nb, "Preflight：执行前校验")
        with pytest.raises(ValueError, match="QC_STRATEGY"):
            exec(src, params_env)

    def test_doublet_threshold_order_inverted(self, tmp_path):
        """DOUBLET_SCORE_LOW > DOUBLET_SCORE_HIGH → ValueError。"""
        import yaml

        manifest_path = tmp_path / "manifest.yaml"
        dummy_h5ad = tmp_path / "dummy.h5ad"
        dummy_h5ad.write_bytes(b"")
        with open(manifest_path, "w") as f:
            yaml.dump({"input": {"format": "10x_h5", "path": str(dummy_h5ad)}, "source_dataset": "test"}, f)

        nb = _load_notebook()
        params_env = _cell(nb, "=== PARAMS ===")
        params_env["MANIFEST_PATH"] = str(manifest_path)
        params_env["DOUBLET_SCORE_LOW"] = 0.8
        params_env["DOUBLET_SCORE_HIGH"] = 0.3
        params_env["yaml"] = yaml
        params_env["Path"] = Path

        src = _source(nb, "Preflight：执行前校验")
        with pytest.raises(ValueError, match="三态阈值区间倒置"):
            exec(src, params_env)

    def test_soupx_enabled_without_raw_path(self, tmp_path):
        """SOUPX_ENABLED=True 且 manifest 无 raw_matrix_path → ValueError。"""
        import yaml

        manifest_path = tmp_path / "manifest.yaml"
        dummy_h5ad = tmp_path / "dummy.h5ad"
        dummy_h5ad.write_bytes(b"")
        with open(manifest_path, "w") as f:
            yaml.dump({"input": {"format": "10x_h5", "path": str(dummy_h5ad)}, "source_dataset": "test"}, f)

        nb = _load_notebook()
        params_env = _cell(nb, "=== PARAMS ===")
        params_env["MANIFEST_PATH"] = str(manifest_path)
        params_env["SOUPX_ENABLED"] = True
        params_env["yaml"] = yaml
        params_env["Path"] = Path

        src = _source(nb, "Preflight：执行前校验")
        with pytest.raises(ValueError, match="raw_matrix_path"):
            exec(src, params_env)


# ---------------------------------------------------------------------------
# G. 科学参数四要素注释在场（静态、软校验）
# ---------------------------------------------------------------------------

class TestFourElementComments:
    """验证 PARAMS cell source 含四要素注释关键词。"""

    def test_four_element_keywords_present(self):
        """PARAMS cell source 含「何时」「默认」「风险」等关键词至少覆盖主要参数。"""
        nb = _load_notebook()
        src = _source(nb, "=== PARAMS ===")
        # 至少覆盖 QC_STRATEGY / N_MAD / SOUPX_ENABLED 三处
        assert "何时" in src, "missing '何时' (when to change)"
        assert "默认依据" in src, "missing '默认依据' (default basis)"
        assert "调大调小" in src or "调高" in src, "missing impact keywords"
        assert "是什么" in src, "missing '是什么' (what it is)"

    def test_qc_strategy_has_four_elements(self):
        """QC_STRATEGY 邻近区域含四要素关键词。"""
        nb = _load_notebook()
        src = _source(nb, "=== PARAMS ===")
        # 从 QC_STRATEGY 赋值行向前找注释
        lines = src.split("\n")
        qc_line = next(i for i, line in enumerate(lines) if line.strip().startswith("QC_STRATEGY"))
        context = "\n".join(lines[max(0, qc_line - 7) : qc_line + 1])
        assert "是什么" in context
        assert "默认" in context
        assert "何时" in context

    def test_soupx_enabled_has_four_elements(self):
        """SOUPX_ENABLED 邻近区域含四要素关键词。"""
        nb = _load_notebook()
        src = _source(nb, "=== PARAMS ===")
        lines = src.split("\n")
        soupx_line = next(i for i, line in enumerate(lines) if line.strip().startswith("SOUPX_ENABLED"))
        context = "\n".join(lines[max(0, soupx_line - 6) : soupx_line + 1])
        assert "是什么" in context
        assert "默认" in context
        assert "何时" in context


# ---------------------------------------------------------------------------
# H. 红线守护（复用 P1-c 断言口径，证明骨架未破坏既成成果）
# ---------------------------------------------------------------------------

class TestRedlineGuards:
    """验证 P1-e 骨架变更未破坏 P0+P1-c 既成成果。"""

    # -- counts 契约 --

    def test_counts_contract_cell_untouched(self):
        """含「建立 counts 契约」的 cell 仍含 layers["counts"] 与 expression_contract，
        且不含任何 doublet 词。"""
        nb = _load_notebook()
        src = _source(nb, "建立 counts 契约")
        assert 'layers["counts"]' in src
        assert "expression_contract" in src
        for term in ["doublet_score", "doublet_class", "doublet_include", "doublet_contract"]:
            assert term not in src, (
                f"counts contract cell should not contain {term}"
            )

    # -- doublet 三态 --

    def test_doublet_cell_has_three_state(self):
        """含「双细胞鉴定」cell 仍含 doublet_class + doublet_include。"""
        nb = _load_notebook()
        # 使用更精确的 marker，避免匹配到 markdown cell 或收尾 cell
        src = _source(nb, "双细胞鉴定（决策8）：三态定级")
        assert "doublet_class" in src
        assert "doublet_include" in src

    def test_doublet_contract_persisted(self):
        """含 uns["doublet_contract"] 且 "method": 的 cell 仍含必要键。"""
        nb = _load_notebook()
        # 精确匹配：找到含 uns["doublet_contract"] 赋值（method 键）的 code cell
        found = False
        for cell in nb["cells"]:
            if cell["cell_type"] != "code":
                continue
            src = "".join(cell["source"])
            if 'uns["doublet_contract"]' in src and '"method":' in src:
                found = True
                for key in ["method", "per_sample_thresholds", "needs_review", "random_seed"]:
                    assert key in src, f"doublet_contract should contain key: {key}"
                break
        assert found, "doublet_contract assignment cell not found"

    # -- checkpoint --

    def test_checkpoint_has_doublet_postconditions(self):
        """Checkpoint cell 含 doublet_columns_present 等 hard_postconditions。"""
        nb = _load_notebook()
        src = _source(nb, "Checkpoint：写入 per-dataset")
        for key in [
            "doublet_columns_present",
            "doublet_contract_present",
        ]:
            assert key in src, f"checkpoint missing postcondition: {key}"
        assert "NEEDS_REVIEW" in src or "doublet_needs_review" in src, (
            "checkpoint should handle NEEDS_REVIEW"
        )

    # -- SoupX --

    def test_soupx_cell_exists(self):
        """含「环境 RNA 校正（SoupX）」cell 仍在且含 SOUPX_ENABLED 分支。"""
        nb = _load_notebook()
        # markdown cell
        _find_cell_by_marker(nb, "环境 RNA 校正（SoupX）", cell_type="markdown")
        # code cell
        src = _source(nb, "环境 RNA 校正（SoupX）—— subprocess")
        assert "SOUPX_ENABLED" in src
