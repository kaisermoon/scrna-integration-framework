"""
B9 叶子测试：07_downstream notebook 版本字面量收敛验证。

目标：
- D01-D13 的 PARAMS cell 单点定义 UPSTREAM_VERSION / OUTPUT_VERSION，
  其余 code cell 不再含散落版本字面量（checkpoint 路径版本 / 裸版本值 / stage 元数据键）。
- D14 仅需 UPSTREAM_VERSION；OUTPUT_VERSION=1（int 遗留）不动，且不跑 R2/R3。

检测器基于 AST 遍历 ast.Constant(str)，不依赖注释内容或 git SHA，CI shallow-clone 安全。
"""

import ast
import json
import os
import re

import pytest

# ============================================================
# 检测器（AST-based，精确判别）
# ============================================================

def _find_version_constants(cell_source: str) -> list[tuple[str, int, str]]:
    """
    对单个 code cell 的 source 解析 AST，返回散落的 checkpoint 版本字面量列表。

    返回：[(line_no, val, reason), ...]

    三条规则，命中任一即判定为散落版本字面量：
    R1: re.search(r'_v\\d+\\.h5ad$', s)  — checkpoint 路径版本
    R2: re.fullmatch(r'v\\d+', s)         — 裸版本值（uns["version"] = "v1" 的右值）
    R3: re.fullmatch(r'\\d+[a-z0-9]*_.+_v\\d+', s) — stage 元数据键（数字开头才算是stage前缀）
          天然排除 cell_type_final_v1/pseudotime_monocle3_v1 这类字母开头的字段名。
    """
    findings = []
    try:
        tree = ast.parse(cell_source)
    except SyntaxError:
        return findings  # 语法错误的 cell 跳过
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            s = node.value
            # R1: checkpoint path version
            if re.search(r'_v\d+\.h5ad$', s):
                findings.append((node.lineno if hasattr(node, 'lineno') else 0, s, 'R1-path'))
            # R2: bare version value
            elif re.fullmatch(r'v\d+', s):
                findings.append((node.lineno if hasattr(node, 'lineno') else 0, s, 'R2-bare'))
            # R3: stage metadata key (starts with digit, ends with _v\d+)
            elif re.fullmatch(r'\d+[a-z0-9]*_.+_v\d+', s):
                findings.append((node.lineno if hasattr(node, 'lineno') else 0, s, 'R3-stage-key'))
    return findings


def detect(source: str) -> bool:
    """检测 source 是否含散落版本字面量。True=有违规。"""
    return len(_find_version_constants(source)) > 0


# ============================================================
# 正负对照测试（证明检测器非空跑）
# ============================================================

def test_detector_flags_real_violations():
    """负例：散落版本字面量必须被抓到（证明非空跑）。"""
    # R1: checkpoint path version
    assert detect('x = "results/07_deg_v1.h5ad"')
    assert detect('x = "results/10_pseudotime_v99.h5ad"')
    # R2: bare version value
    assert detect('adata.uns["version"] = "v1"')
    assert detect("adata.uns['version'] = 'v2'")
    # R3: stage metadata key
    assert detect('adata.uns["15_gene_modules_v1"] = {}')
    assert detect("adata.uns['08_pseudobulk_deg_v3'] = {}")


def test_detector_ignores_valid_and_field_names():
    """正例：收敛后写法、字段名契约、局部变量，全部放行（证明无误报）。"""
    # 收敛后写法（f-string 常量段不含 _v\d+.h5ad）
    assert not detect('x = f"results/07_deg_{OUTPUT_VERSION}.h5ad"')
    assert not detect('x = f"results/06_annotated_{UPSTREAM_VERSION}.h5ad"')
    # uns["version"] 引用变量
    assert not detect('adata.uns["version"] = OUTPUT_VERSION')
    # uns stage key 引用变量
    assert not detect('adata.uns[f"07_deg_{OUTPUT_VERSION}"] = x')
    # 字段名契约（字母开头，R3 不命中）
    assert not detect('adata.obs["cell_type_final_v1"]')
    assert not detect('_cand = "pseudotime_monocle3_v1"')
    assert not detect('"cellrank2_pseudotime"')
    assert not detect('"dpt_pseudotime"')
    # 局部变量（非 str 常量 "v1" 的 Python 变量名，R2 是 re.fullmatch 只匹配字符串内容）
    # _v1 如果出现为字符串则测 R2，但 D09 c#4 中 _v1 是 Python 变量名，不触发
    assert not detect("_v1 = proportion_df.loc[mask]")
    assert not detect("_v0 = proportion_df.loc[mask]")
    # 方法标签字典键
    assert not detect('"cytotrace_v1"')


# ==========================================================
# Notebook 定位工具
# ==========================================================

NOTEBOOK_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "notebooks", "07_downstream"
)

D01_D14_NBS = [f"D{i:02d}" for i in range(1, 15)]
# 实际文件名映射
NB_FILES = {
    "D01": "D01_deg.ipynb",
    "D02": "D02_pseudobulk_deg.ipynb",
    "D03": "D03_cnv.ipynb",
    "D04": "D04_pseudotime.ipynb",
    "D05": "D05_pseudotime_monocle3.ipynb",
    "D06": "D06_pseudotime_cellrank2.ipynb",
    "D07": "D07_potency_cytotrace2.ipynb",
    "D08": "D08_pseudotime_compare.ipynb",
    "D09": "D09_abundance.ipynb",
    "D10": "D10_pathway.ipynb",
    "D11": "D11_grn.ipynb",
    "D12": "D12_cell_communication.ipynb",
    "D13": "D13_gene_modules.ipynb",
    "D14": "D14_trajectory_de.ipynb",
}


def _load_nb(nb_id: str) -> dict:
    """加载 notebook JSON。"""
    path = os.path.join(NOTEBOOK_DIR, NB_FILES[nb_id])
    with open(path) as f:
        return json.load(f)


def _code_cells(nb: dict) -> list[tuple[int, str]]:
    """
    返回所有 code cell 的 (code_cell_index, source) 列表。
    code_cell_index 是 1-based（只计 code cell，不计 markdown）。
    """
    result = []
    idx = 0
    for cell in nb['cells']:
        if cell['cell_type'] == 'code':
            idx += 1
            result.append((idx, ''.join(cell['source'])))
    return result


def _has_ast_assign_str(source: str, var_name: str) -> bool:
    """
    用 AST 验证 source 中存在对 var_name 的真赋值（value 为 str 常量）。
    排除在注释/字符串内误匹配（P3-c 教训）。
    """
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return False
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if (isinstance(target, ast.Name) and target.id == var_name
                        and isinstance(node.value, ast.Constant)
                        and isinstance(node.value.value, str)):
                    return True
    return False


# ==========================================================
# 主断言：D01–D13
# ==========================================================

@pytest.mark.parametrize("nb_id", [n for n in D01_D14_NBS if n != "D14"])
def test_params_cell_defines_version_vars(nb_id: str):
    """
    PARAMS cell（第一个 code cell）用 AST Assign 定义了 UPSTREAM_VERSION
    和 OUTPUT_VERSION，value 均为 str 常量。

    为什么用 AST 而非字符串搜索：避免注释中的变量名产生假阳性（P3-c 教训）。
    """
    nb = _load_nb(nb_id)
    _, params_source = _code_cells(nb)[0]  # first code cell

    assert _has_ast_assign_str(params_source, 'UPSTREAM_VERSION'), \
        f"{nb_id}: PARAMS cell 中缺少 str 赋值 `UPSTREAM_VERSION = \"v1\"`"
    assert _has_ast_assign_str(params_source, 'OUTPUT_VERSION'), \
        f"{nb_id}: PARAMS cell 中缺少 str 赋值 `OUTPUT_VERSION = \"v1\"`"


@pytest.mark.parametrize("nb_id", [n for n in D01_D14_NBS if n != "D14"])
def test_no_scattered_version_literals_outside_params(nb_id: str):
    """
    除 PARAMS cell 之外的所有 code cell 均不含散落版本字面量。
    检测器基于 AST Constant(str)，R1+R2+R3 三条规则。

    为什么排除 PARAMS：PARAMS 中的 "v1" 是版本变量的默认值本身，
    允许在那里出现。其余所有 cell 必须通过检测器。
    """
    nb = _load_nb(nb_id)
    cells = _code_cells(nb)

    violations = []
    for ci, source in cells[1:]:  # skip first code cell (PARAMS)
        findings = _find_version_constants(source)
        for lineno, val, reason in findings:
            violations.append(f"  c#[{ci}] L{lineno}: {val!r} ({reason})")

    assert not violations, \
        f"{nb_id}: 非 PARAMS cell {len(violations)} 处散落版本字面量:\n" + "\n".join(violations)


# ==========================================================
# D14 特例分支
# ==========================================================

def test_d14_params_defines_upstream_version():
    """
    D14 的 PARAMS cell 定义了 UPSTREAM_VERSION（str 常量）。
    不强制要求 OUTPUT_VERSION（遗留 int=1，不动）。
    """
    nb = _load_nb("D14")
    _, params_source = _code_cells(nb)[0]

    assert _has_ast_assign_str(params_source, 'UPSTREAM_VERSION'), \
        "D14: PARAMS cell 中缺少 str 赋值 `UPSTREAM_VERSION = \"v1\"`"


def test_d14_no_r1_outside_params():
    """
    D14 非 PARAMS cell 无 checkpoint 路径版本字面量（R1）。
    不跑 R2/R3（D14 不产出带版本号的 checkpoint h5ad，无 uns["version"] 写入）。
    """
    nb = _load_nb("D14")
    cells = _code_cells(nb)

    violations = []
    for ci, source in cells[1:]:
        for lineno, val, reason in _find_version_constants(source):
            if reason == 'R1-path':
                violations.append(f"  c#[{ci}] L{lineno}: {val!r} ({reason})")

    assert not violations, \
        f"D14: 非 PARAMS cell {len(violations)} 处 checkpoint 路径版本字面量:\n" + "\n".join(violations)
