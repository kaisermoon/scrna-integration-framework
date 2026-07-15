"""测试 P1-e-nancang：UX 骨架叶子（PARAMS 四组化 + preflight 校验 + 科学参数四要素注释）。

所有测试静态检查 notebooks/01_per_dataset/01_nancang.ipynb 的 cell 结构与源码，
T6 通过 exec 验证 preflight 校验行为（不依赖 nbconvert 端到端）。
"""
from __future__ import annotations

import json
import os
import re
import sys
import tempfile
from pathlib import Path

import pytest

# ---- helpers ----------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
_NB_PATH = ROOT / "notebooks" / "01_per_dataset" / "01_nancang.ipynb"


def _nb() -> dict:
    """加载 notebook JSON。"""
    assert _NB_PATH.is_file(), f"notebook not found: {_NB_PATH}"
    with open(_NB_PATH, encoding="utf-8") as f:
        return json.load(f)


def _source_cells(marker: str, cell_type: str = "code") -> list[tuple[int, str]]:
    """返回所有包含 marker 的 cell 的 (index, source) 列表。"""
    results = []
    for i, cell in enumerate(_nb()["cells"]):
        if cell["cell_type"] != cell_type:
            continue
        src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
        if marker in src:
            results.append((i, src))
    return results


def _first_source(marker: str, cell_type: str = "code") -> str:
    """返回首个包含 marker 的 cell 的 source 字符串。"""
    results = _source_cells(marker, cell_type)
    assert results, f"No {cell_type} cell with marker {marker!r} found"
    return results[0][1]


def _cell_index(marker: str, cell_type: str = "code") -> int:
    """返回首个包含 marker 的 cell 的 index。"""
    results = _source_cells(marker, cell_type)
    assert results, f"No {cell_type} cell with marker {marker!r} found"
    return results[0][0]


def _cell_md_index(marker: str) -> int:
    """返回首个包含 marker 的 markdown cell 的 index。"""
    return _cell_index(marker, "markdown")


def _collect_code_sources() -> dict[int, str]:
    """返回所有 code cell 的 {index: source} 映射。"""
    result = {}
    for i, cell in enumerate(_nb()["cells"]):
        if cell["cell_type"] != "code":
            continue
        src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
        result[i] = src
    return result


# ---- T1: PARAMS 四组 markdown 小标题存在且顺序正确 --------------------------

def test_t1_group_headers_order():
    """四组 markdown 小标题存在且顺序为 数据源 → QC 阈值 → 方法开关 → 输出版本与运行标识。"""
    nb = _nb()
    expected = [
        "### 1. 数据源",
        "### 2. QC 阈值",
        "### 3. 方法开关",
        "### 4. 输出版本与运行标识",
    ]
    found = []
    for cell in nb["cells"]:
        if cell["cell_type"] != "markdown":
            continue
        src = cell["source"] if isinstance(cell["source"], str) else "".join(cell["source"])
        for exp in expected:
            if exp in src:
                found.append(exp)
    assert found == expected, f"Group headers order mismatch. Expected {expected}, got {found}"


# ---- T2: 默认值语义未变 -----------------------------------------------------

def test_t2_default_values_preserved():
    """解析 PARAMS 四个 code cell 源码，断言关键字面量原样。"""
    # 收集四个 group 的 code cell
    g1_idx = _cell_md_index("### 1. 数据源")
    g2_idx = _cell_md_index("### 2. QC 阈值")
    g3_idx = _cell_md_index("### 3. 方法开关")
    g4_idx = _cell_md_index("### 4. 输出版本与运行标识")

    nb = _nb()
    # code cells are immediately after their md headers
    g1_src = "".join(nb["cells"][g1_idx + 1]["source"])
    g2_src = "".join(nb["cells"][g2_idx + 1]["source"])
    g3_src = "".join(nb["cells"][g3_idx + 1]["source"])
    g4_src = "".join(nb["cells"][g4_idx + 1]["source"])
    all_src = g1_src + "\n" + g2_src + "\n" + g3_src + "\n" + g4_src

    # 逐字断言关键默认值
    checks = [
        ('SOUPX_ENABLED = True', "SOUPX_ENABLED 默认值"),
        ('SOUPX_INTEGER_ROUND = True', "SOUPX_INTEGER_ROUND 默认值"),
        ('N_MAD = 5', "N_MAD 默认值"),
        ('EXPECTED_DOUBLET_RATE = 0.06', "EXPECTED_DOUBLET_RATE 默认值"),
        ('DOUBLET_RATE_ALERT_HIGH = 0.40', "DOUBLET_RATE_ALERT_HIGH 默认值"),
        ('DOUBLET_MIN_CELLS = 50', "DOUBLET_MIN_CELLS 默认值"),
        ('MIN_CELLS_PER_GENE = 3', "MIN_CELLS_PER_GENE 默认值"),
        ('QC_STRATEGY = "adaptive"', "QC_STRATEGY 默认值"),
        ('OUTPUT_VERSION = 1', "OUTPUT_VERSION 默认值"),
        ('RANDOM_SEED = 42', "RANDOM_SEED 默认值"),
        ('DOUBLET_UNCERTAIN_INCLUDE = True', "DOUBLET_UNCERTAIN_INCLUDE 默认值"),
    ]
    for expected, desc in checks:
        assert expected in all_src, f"{desc} 未找到: '{expected}'"

    # 校验所有原 PARAMS 变量名均被赋值（不丢失）
    required_vars = [
        "MANIFEST_PATH", "RUN_ID", "RUN_ROOT", "OUTPUT_FILENAME",
        "QC_STRATEGY", "N_MAD", "RANDOM_SEED",
        "SOUPX_ENABLED", "SOUPX_INTEGER_ROUND", "EXPECTED_DOUBLET_RATE",
        "FLAG_HEMOGLOBIN", "HB_THRESHOLD_PCT", "FLAG_STRESS_GENES",
        "SCORE_CELL_CYCLE", "DOUBLET_SCORE_THRESHOLD", "DOUBLET_SCORE_HIGH",
        "DOUBLET_SCORE_LOW", "DOUBLET_UNCERTAIN_INCLUDE",
        "DOUBLET_RATE_ALERT_HIGH", "DOUBLET_RATE_ALERT_LOW", "DOUBLET_MIN_CELLS",
        "OUTPUT_VERSION", "PER_SAMPLE_MAD", "MIN_CELLS_PER_GENE",
        "MIN_GENES", "MAX_GENES", "MIN_COUNTS", "MAX_PCT_MT",
    ]
    for var in required_vars:
        assert re.search(rf'\b{var}\s*=', all_src), f"变量 {var} 未被赋值"


# ---- T3: preflight cell 位置与内容正确 -------------------------------------

def test_t3_preflight_position():
    """preflight cell 位于 Setup 之后、数据读入之前。"""
    setup_idx = _cell_index("# === Setup")
    data_idx = _cell_index("# 数据读入：10x mtx")

    nb = _nb()

    # preflight md 应在 setup+1 或 setup+2 位置
    pf_md_idx = _cell_md_index("Preflight：运行前配置校验")
    pf_code_idx = _cell_index("# Preflight：运行前配置校验")

    assert setup_idx < pf_md_idx < data_idx, (
        f"preflight md ({pf_md_idx}) 不在 Setup ({setup_idx}) 和 data loading ({data_idx}) 之间"
    )
    assert setup_idx < pf_code_idx < data_idx, (
        f"preflight code ({pf_code_idx}) 不在 Setup ({setup_idx}) 和 data loading ({data_idx}) 之间"
    )
    assert pf_md_idx + 1 == pf_code_idx, "preflight md 和 code 应相邻"


def test_t3_preflight_contains_checks():
    """preflight code cell 含对 MANIFEST_PATH / SOUPX_ENABLED / soupx_run.R 校验及至少一处 raise。"""
    pf_src = _first_source("# Preflight：运行前配置校验")

    assert "MANIFEST_PATH" in pf_src, "preflight 应校验 MANIFEST_PATH"
    assert "SOUPX_ENABLED" in pf_src, "preflight 应校验 SOUPX_ENABLED"
    assert "soupx_run.R" in pf_src, "preflight 应校验 soupx_run.R"
    assert "raise" in pf_src, "preflight 应含至少一处 raise（校验失败抛异常）"


# ---- T4: redline marker 完好 --------------------------------------------------

def test_t4_redline_count_contract():
    """建立 layers['counts'] 的行存在且逻辑未变。"""
    src = _first_source('layers["counts"] = sp.csr_matrix(adata.X')
    assert "sp.csr_matrix(adata.X" in src
    assert 'dtype=np.float32' in src


def test_t4_redline_counts_checksum():
    """_counts_checksum 定义行存在。"""
    src = _first_source("_counts_checksum")
    assert "_counts_checksum = float" in src
    assert "_counts_checksum_nnz = int" in src


def test_t4_redline_counts_soupx_write():
    """layers['counts_soupx'] 写入存在。"""
    result = _source_cells('layers["counts_soupx"]')
    assert len(result) >= 1, "counts_soupx 写入 cell 必须存在"


def test_t4_redline_doublet_class():
    """doublet_class 三态列写入存在。"""
    result = _source_cells("doublet_class")
    assert len(result) >= 1, "doublet_class 写入 cell 必须存在"


def test_t4_redline_doublet_include():
    """doublet_include 派生存在。"""
    result = _source_cells("doublet_include")
    assert len(result) >= 1, "doublet_include 派生 cell 必须存在"


def test_t4_redline_doublet_contract():
    """doublet_contract 写入存在。"""
    result = _source_cells("doublet_contract")
    assert len(result) >= 1, "doublet_contract cell 必须存在"


def test_t4_redline_checkpoint_hard_postconditions():
    """checkpoint 的 hard_postconditions 存在。"""
    result = _source_cells("hard_postconditions")
    assert len(result) >= 1, "hard_postconditions cell 必须存在"


def test_t4_redline_doublet_hd_guard():
    """checkpoint 的 _hd guard 存在。"""
    src = _first_source("hard_postconditions")
    assert '_hd = "doublet_class" in adata.obs.columns' in src


def test_t4_redline_counts_layer_unchanged():
    """counts_layer_unchanged guard 存在。"""
    src = _first_source("hard_postconditions")
    assert "counts_layer_unchanged" in src


def test_t4_redline_expression_contract():
    """expression_contract 8字段 schema 存在。"""
    result = _source_cells("expression_contract")
    assert len(result) >= 1, "expression_contract cell 必须存在"


# ---- T5: 科学参数四要素注释 --------------------------------------------------

def test_t5_four_element_comments():
    """组 2（QC 阈值）和组 3（方法开关）的科学参数上方含四要素注释信号。"""
    g2_idx = _cell_md_index("### 2. QC 阈值")
    g3_idx = _cell_md_index("### 3. 方法开关")

    nb = _nb()
    g2_src = "".join(nb["cells"][g2_idx + 1]["source"])
    g3_src = "".join(nb["cells"][g3_idx + 1]["source"])

    # 科学参数列表
    sci_params_g2 = [
        "QC_STRATEGY", "N_MAD", "PER_SAMPLE_MAD", "MIN_CELLS_PER_GENE",
        "MIN_GENES", "MAX_GENES", "MIN_COUNTS", "MAX_PCT_MT", "HB_THRESHOLD_PCT",
    ]
    sci_params_g3 = [
        "SOUPX_ENABLED", "SOUPX_INTEGER_ROUND", "EXPECTED_DOUBLET_RATE",
        "DOUBLET_SCORE_THRESHOLD", "DOUBLET_SCORE_HIGH", "DOUBLET_SCORE_LOW",
        "DOUBLET_UNCERTAIN_INCLUDE", "DOUBLET_RATE_ALERT_HIGH",
        "DOUBLET_RATE_ALERT_LOW", "DOUBLET_MIN_CELLS",
        "FLAG_HEMOGLOBIN", "FLAG_STRESS_GENES", "SCORE_CELL_CYCLE",
    ]

    four_element_signals = ["是什么", "默认依据", "默认", "调", "何时该改", "何时"]

    for param in sci_params_g2:
        _check_param_has_block_comment(param, g2_src, four_element_signals)

    for param in sci_params_g3:
        _check_param_has_block_comment(param, g3_src, four_element_signals)


def _check_param_has_block_comment(param: str, src: str, signals: list[str]):
    """校验参数赋值语句上方有块注释（含四要素信号中的至少 2 个）。"""
    # 找到参数赋值位置
    match = re.search(rf'^{param}\s*=\s*', src, re.MULTILINE)
    assert match is not None, f"参数 {param} 未在对应组中找到赋值"
    pos = match.start()

    # 取赋值行之前的内容（注释块应在前面）
    before = src[:pos]
    # 取最后 15 行作为可能的注释块
    before_lines = before.split("\n")
    comment_block = "\n".join(before_lines[-15:])

    # 检查是否含至少 2 个四要素信号词
    hit_count = sum(1 for sig in signals if sig in comment_block)
    assert hit_count >= 2, (
        f"参数 {param} 的四要素注释不足（命中信号 {hit_count}/2+）。"
        f"注释块末尾:\n{comment_block[-400:]}"
    )


# ---- T6: preflight 功能测试（真的 raise，不是摆设）--------------------------

class TestPreflightFunctional:
    """通过 exec 验证 preflight 校验在非法配置时 raise、正常配置不 raise。"""

    @staticmethod
    def _extract_group_code() -> dict[str, str]:
        """提取四个分组 code cell 的源码。"""
        nb = _nb()
        g1_idx = _cell_md_index("### 1. 数据源")
        g2_idx = _cell_md_index("### 2. QC 阈值")
        g3_idx = _cell_md_index("### 3. 方法开关")
        g4_idx = _cell_md_index("### 4. 输出版本与运行标识")
        return {
            "g1": "".join(nb["cells"][g1_idx + 1]["source"]),
            "g2": "".join(nb["cells"][g2_idx + 1]["source"]),
            "g3": "".join(nb["cells"][g3_idx + 1]["source"]),
            "g4": "".join(nb["cells"][g4_idx + 1]["source"]),
        }

    @staticmethod
    def _extract_preflight_code() -> str:
        """提取 preflight code cell 源码。"""
        return _first_source("# Preflight：运行前配置校验")

    @staticmethod
    def _make_exec_namespace(tmp_path: Path, **overrides):
        """构建 preflight exec 用的命名空间。"""
        groups = TestPreflightFunctional._extract_group_code()
        pf_src = TestPreflightFunctional._extract_preflight_code()

        # 拼合：先 exec 四组 PARAMS，再 exec preflight
        combined = (
            groups["g1"] + "\n" + groups["g2"] + "\n" +
            groups["g3"] + "\n" + groups["g4"] + "\n"
        )
        ns: dict = {}
        exec(combined, ns)

        # 应用 overrides（如 N_MAD=-1 测试非法值）
        for k, v in overrides.items():
            ns[k] = v

        # 注入 preflight 依赖的外部变量（由 Setup cell 提供）
        ns["RSCRIPT_BIN"] = "Rscript"
        ns["R_AVAILABLE"] = True

        return ns, pf_src

    def test_t6a_normal_config(self, tmp_path):
        """正常配置（manifest.yaml + input dir + raw dir + soupx_run.R 均存在）→ 不 raise。"""
        # 造 manifest
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        (data_dir / "matrix.mtx.gz").touch()

        raw_dir = tmp_path / "raw_dir"
        raw_dir.mkdir()

        manifest = {
            "input": {"path": str(data_dir), "format": "10x_mtx", "raw_path": str(raw_dir)},
            "source_dataset": "nancang",
        }
        import yaml
        manifest_path = tmp_path / "manifest.yaml"
        with open(manifest_path, "w") as f:
            yaml.dump(manifest, f)

        # 造 soupx_run.R
        scripts_dir = tmp_path / "scripts"
        scripts_dir.mkdir()
        (scripts_dir / "soupx_run.R").touch()

        # 构建命名空间 + 切工作目录
        ns, pf_src = self._make_exec_namespace(
            tmp_path, MANIFEST_PATH=str(manifest_path),
        )

        import os as _os
        _old_cwd = _os.getcwd()
        try:
            _os.chdir(str(tmp_path))
            exec(pf_src, ns)
        finally:
            _os.chdir(_old_cwd)

    def test_t6b_manifest_not_found(self, tmp_path):
        """MANIFEST_PATH 指向不存在文件 → raise。"""
        ns, pf_src = self._make_exec_namespace(
            tmp_path, MANIFEST_PATH=str(tmp_path / "nonexistent.yaml"),
        )
        with pytest.raises((FileNotFoundError, SystemExit)):
            exec(pf_src, ns)

    def test_t6c_n_mad_negative(self, tmp_path):
        """N_MAD=-1 → raise AssertionError。"""
        # 先造正常 manifest 环境
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        (data_dir / "matrix.mtx.gz").touch()

        manifest = {
            "input": {"path": str(data_dir), "format": "10x_mtx"},
            "source_dataset": "nancang",
        }
        import yaml
        manifest_path = tmp_path / "manifest.yaml"
        with open(manifest_path, "w") as f:
            yaml.dump(manifest, f)

        ns, pf_src = self._make_exec_namespace(
            tmp_path, MANIFEST_PATH=str(manifest_path), N_MAD=-1, SOUPX_ENABLED=False,
        )
        import os as _os
        _old_cwd = _os.getcwd()
        try:
            _os.chdir(str(tmp_path))
            with pytest.raises(AssertionError):
                exec(pf_src, ns)
        finally:
            _os.chdir(_old_cwd)

    def test_t6d_soupx_no_raw_dir(self, tmp_path):
        """SOUPX_ENABLED=True 但 raw_path 目录不存在 → raise。"""
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        (data_dir / "matrix.mtx.gz").touch()

        raw_dir = tmp_path / "raw_dir_nonexistent"  # 不创建！

        manifest = {
            "input": {"path": str(data_dir), "format": "10x_mtx", "raw_path": str(raw_dir)},
            "source_dataset": "nancang",
        }
        import yaml
        manifest_path = tmp_path / "manifest.yaml"
        with open(manifest_path, "w") as f:
            yaml.dump(manifest, f)

        scripts_dir = tmp_path / "scripts"
        scripts_dir.mkdir()
        (scripts_dir / "soupx_run.R").touch()

        ns, pf_src = self._make_exec_namespace(tmp_path, MANIFEST_PATH=str(manifest_path))
        import os as _os
        _old_cwd = _os.getcwd()
        try:
            _os.chdir(str(tmp_path))
            with pytest.raises(FileNotFoundError):
                exec(pf_src, ns)
        finally:
            _os.chdir(_old_cwd)

    def test_t6e_soupx_no_r_script(self, tmp_path):
        """SOUPX_ENABLED=True 但 soupx_run.R 不存在 → raise。"""
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        (data_dir / "matrix.mtx.gz").touch()

        raw_dir = tmp_path / "raw_dir"
        raw_dir.mkdir()

        manifest = {
            "input": {"path": str(data_dir), "format": "10x_mtx", "raw_path": str(raw_dir)},
            "source_dataset": "nancang",
        }
        import yaml
        manifest_path = tmp_path / "manifest.yaml"
        with open(manifest_path, "w") as f:
            yaml.dump(manifest, f)

        # 不创建 scripts/soupx_run.R

        ns, pf_src = self._make_exec_namespace(tmp_path, MANIFEST_PATH=str(manifest_path))
        import os as _os
        _old_cwd = _os.getcwd()
        try:
            _os.chdir(str(tmp_path))
            with pytest.raises(FileNotFoundError):
                exec(pf_src, ns)
        finally:
            _os.chdir(_old_cwd)

    def test_t6f_soupx_r_not_available(self, tmp_path):
        """SOUPX_ENABLED=True 且 R_AVAILABLE=False → 不 raise（只 WARNING）。"""
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        (data_dir / "matrix.mtx.gz").touch()

        raw_dir = tmp_path / "raw_dir"
        raw_dir.mkdir()

        manifest = {
            "input": {"path": str(data_dir), "format": "10x_mtx", "raw_path": str(raw_dir)},
            "source_dataset": "nancang",
        }
        import yaml
        manifest_path = tmp_path / "manifest.yaml"
        with open(manifest_path, "w") as f:
            yaml.dump(manifest, f)

        scripts_dir = tmp_path / "scripts"
        scripts_dir.mkdir()
        (scripts_dir / "soupx_run.R").touch()

        ns, pf_src = self._make_exec_namespace(tmp_path, MANIFEST_PATH=str(manifest_path))
        ns["R_AVAILABLE"] = False  # 模拟 R 不可用

        import os as _os
        _old_cwd = _os.getcwd()
        try:
            _os.chdir(str(tmp_path))
            exec(pf_src, ns)  # 不应 raise
        finally:
            _os.chdir(_old_cwd)
