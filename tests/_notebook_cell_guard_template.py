"""
notebook cell 守卫测试工具函数
==============================
这些工具函数用来防止 notebook cell 守卫测试的三类空跑陷阱：

1. AST 搜索失配：用 cell id 定位 cell 而 id 随版本漂移，导致找不到目标 cell 但测试静默通过
2. 断言目标错位：断言对象不是真正的 gate 变量，测试通过但没有验证任何有意义的约束
3. 同语反复：断言形如 list(x) == list(x)，永远为真，测试毫无判别力

用法：在 test_*.py 中直接 import 并调用，不是 pytest fixture。
"""

import ast
import json
import re
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# 底层 notebook 读取
# ---------------------------------------------------------------------------


def _load_notebook(nb_path: str | Path) -> dict:
    """
    加载 Jupyter notebook 文件，返回解析后的 dict。

    参数：
        nb_path: notebook 文件的绝对路径（.ipynb）
    """
    path = Path(nb_path)
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def _cell_source(cell: dict) -> str:
    """
    将 cell 的 source 字段（list[str] 或 str）合并为单个字符串。

    参数：
        cell: notebook cell dict
    """
    src = cell.get("source", "")
    if isinstance(src, list):
        return "".join(src)
    return src


# ---------------------------------------------------------------------------
# 公共工具函数
# ---------------------------------------------------------------------------


def find_cell_by_source_fragment(nb_path: str | Path, fragment: str) -> dict | None:
    """
    在 notebook 中查找第一个 source 包含指定子串的 cell。

    比用 cell id 定位更稳定，因为代码内容比 cell id 更少随版本漂移。

    参数：
        nb_path:  notebook 文件绝对路径
        fragment: 要搜索的代码片段字符串（子串匹配）

    返回：
        第一个匹配的 cell dict；找不到返回 None。
    """
    nb = _load_notebook(nb_path)
    for cell in nb.get("cells", []):
        if fragment in _cell_source(cell):
            return cell
    return None


def find_gate_cell(nb_path: str | Path, gate_var_name: str) -> dict | None:
    """
    在 notebook 中查找包含 gate_var_name 赋值语句的 cell。

    匹配模式为 `gate_var_name =` 或 `gate_var_name=`（允许赋值号两侧有空格）。

    参数：
        nb_path:       notebook 文件绝对路径
        gate_var_name: gate 变量名，例如 "PI_CONFIRMED"

    返回：
        第一个匹配的 cell dict；找不到返回 None。
    """
    # 匹配赋值语句：变量名后跟可选空格再跟 `=`（不匹配 `==`）
    pattern = re.compile(r"\b" + re.escape(gate_var_name) + r"\s*=[^=]")
    nb = _load_notebook(nb_path)
    for cell in nb.get("cells", []):
        if pattern.search(_cell_source(cell)):
            return cell
    return None


def extract_gate_value_from_cell(cell: dict, gate_var_name: str) -> Any:
    """
    从 cell source 中解析 gate_var_name 的赋值值。

    支持三种模式：
    - ``gate_var_name = True``
    - ``gate_var_name = False``
    - ``gate_var_name = not <expr>``

    参数：
        cell:          notebook cell dict
        gate_var_name: gate 变量名

    返回：
        True / False / "not <expr>"（字符串形式）；解析失败返回 None。
    """
    source = _cell_source(cell)
    # 简单行级扫描，逐行尝试解析赋值
    pattern = re.compile(
        r"^\s*" + re.escape(gate_var_name) + r"\s*=\s*(.+)$",
        re.MULTILINE,
    )
    match = pattern.search(source)
    if not match:
        return None

    rhs = match.group(1).strip()

    if rhs == "True":
        return True
    if rhs == "False":
        return False
    if rhs.startswith("not "):
        return rhs  # 返回字符串形式，如 "not some_condition"

    # 尝试 AST 字面量求值（处理括号等变体）
    try:
        return ast.literal_eval(rhs)
    except (ValueError, SyntaxError):
        return None


def assert_cell_exists(
    nb_path: str | Path, fragment: str, msg: str | None = None
) -> None:
    """
    断言 notebook 中存在包含 fragment 子串的 cell。

    参数：
        nb_path:  notebook 文件绝对路径
        fragment: 必须出现的代码片段
        msg:      断言失败时的自定义消息（可选）

    异常：
        AssertionError：找不到包含 fragment 的 cell 时抛出。
    """
    cell = find_cell_by_source_fragment(nb_path, fragment)
    if cell is None:
        default_msg = (
            f"notebook {nb_path!r} 中未找到包含以下片段的 cell：\n{fragment!r}"
        )
        raise AssertionError(msg or default_msg)


def assert_gate_not_hardcoded_true(
    nb_path: str | Path, gate_var_name: str
) -> None:
    """
    断言 gate 变量不是简单硬编码为 True。

    防止 coder 在调试后忘记还原，留下 ``PI_CONFIRMED = True`` 使安全闸失效。

    参数：
        nb_path:       notebook 文件绝对路径
        gate_var_name: gate 变量名，例如 "PI_CONFIRMED"

    异常：
        AssertionError：gate cell 不存在，或 gate 变量被硬编码为 True 时抛出。
    """
    cell = find_gate_cell(nb_path, gate_var_name)
    if cell is None:
        raise AssertionError(
            f"notebook {nb_path!r} 中未找到 {gate_var_name!r} 的赋值 cell。"
            " 请检查变量名是否正确。"
        )

    value = extract_gate_value_from_cell(cell, gate_var_name)
    if value is True:
        raise AssertionError(
            f"{gate_var_name!r} 在 notebook {nb_path!r} 中被硬编码为 True，"
            " 安全闸已失效。请将其改为条件表达式或 False。"
        )


# ---------------------------------------------------------------------------
# 同语反复检测
# ---------------------------------------------------------------------------


def assert_no_tautology(test_file_path: str | Path, test_func_name: str) -> None:
    """
    扫描测试函数源码，检测典型的同语反复断言模式。

    检测以下模式（简单 AST 扫描，不覆盖所有情形）：
    - ``list(x) == list(x)``：相同表达式两侧转 list 比较
    - ``x is x``：对象与自身做 identity 检查
    - ``len(x) == len(x)``：相同表达式取 len 后比较

    参数：
        test_file_path: 测试文件的绝对路径
        test_func_name: 要检查的测试函数名（不含括号）

    异常：
        AssertionError：发现同语反复模式时抛出，消息中说明具体模式与位置。
        ValueError：找不到指定函数时抛出。
    """
    source = Path(test_file_path).read_text(encoding="utf-8")
    tree = ast.parse(source)

    # 找到目标函数定义节点
    func_node = None
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name == test_func_name:
                func_node = node
                break

    if func_node is None:
        raise ValueError(
            f"在 {test_file_path!r} 中未找到函数 {test_func_name!r}。"
        )

    tautologies_found = []

    for node in ast.walk(func_node):
        # 模式 1：x is x（同一节点做 is 比较）
        if isinstance(node, ast.Compare):
            left = node.left
            for op, comparator in zip(node.ops, node.comparators, strict=True):
                if isinstance(op, ast.Is):
                    if ast.dump(left) == ast.dump(comparator):
                        tautologies_found.append(
                            f"第 {node.lineno} 行：`x is x` 同语反复"
                            f"（左右表达式相同：{ast.unparse(left)!r}）"
                        )

                # 模式 2：list(x) == list(x) 或 len(x) == len(x)
                if isinstance(op, (ast.Eq, ast.NotEq)):
                    left_dump = ast.dump(left)
                    right_dump = ast.dump(comparator)
                    if left_dump == right_dump:
                        # 是否是包装调用形式
                        if isinstance(left, ast.Call):
                            func_name = _get_call_func_name(left)
                            if func_name in ("list", "len", "set", "tuple", "sorted"):
                                tautologies_found.append(
                                    f"第 {node.lineno} 行："
                                    f"`{func_name}(x) == {func_name}(x)` 同语反复"
                                    f"（表达式：{ast.unparse(node)!r}）"
                                )
                        else:
                            # 纯同名变量比较 x == x
                            if isinstance(left, ast.Name):
                                tautologies_found.append(
                                    f"第 {node.lineno} 行："
                                    f"`{left.id} == {left.id}` 同语反复"
                                )

    if tautologies_found:
        detail = "\n".join(f"  - {t}" for t in tautologies_found)
        raise AssertionError(
            f"函数 {test_func_name!r} 中发现同语反复断言：\n{detail}"
        )


def _get_call_func_name(call_node: ast.Call) -> str | None:
    """
    从 ast.Call 节点中提取被调用函数的名称（仅处理简单名称，不处理属性调用）。

    参数：
        call_node: ast.Call 节点

    返回：
        函数名字符串；无法提取时返回 None。
    """
    if isinstance(call_node.func, ast.Name):
        return call_node.func.id
    return None
