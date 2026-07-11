"""启动脚手架：项目根查找 + 环境变量设置 + sys.path 注入。

本模块是技术管道，不含分析逻辑。每个 notebook 的启动 cell 一行调用即可，
不必重复写 20 行样板代码。遵循 "platform 决策可见" 原则——verbose 模式
print 出所有决策结果（项目根路径、线程数、sys.path 添加情况），
让 notebook cell 输出区直接显示环境诊断。

使用方式（notebook 启动 cell）::

    from scrna_integration.bootstrap import init
    _root = init()
    # 一行替代原来 20 行样板代码

面向非计算机专业学生：注释解释每一步在做什么、为什么需要这样做。
"""

from __future__ import annotations

import os
import sys


def init(verbose: bool = True) -> str:
    """初始化项目运行环境，返回项目根绝对路径。

    按顺序执行三件事：
    1. 找到项目根（含 ``src/scrna_integration/`` 的目录）
    2. 设置 BLAS/线程相关环境变量（防止服务器核心数爆炸）
    3. 将 ``{项目根}/src`` 加入 ``sys.path``（使 ``import scrna_integration`` 可用）

    Parameters
    ----------
    verbose : bool
        True 时 print 出所有决策结果（项目根路径、线程数、sys.path）。
        遵循 "platform 决策可见" 原则——notebook cell 输出区直接显示。

    Returns
    -------
    str
        项目根绝对路径。
    """
    # ---- 1. 找项目根 ----
    # 从当前工作目录向上逐级查找，直到找到包含 src/scrna_integration 的目录。
    # 兼容任意嵌套深度（notebook 可能在 notebooks/01_per_dataset/ 下，
    # 也可能被拷贝到其他地方运行）。
    _root = os.getcwd()
    while _root != os.path.dirname(_root):
        if os.path.isdir(os.path.join(_root, "src", "scrna_integration")):
            break
        _root = os.path.dirname(_root)

    if not os.path.isdir(os.path.join(_root, "src", "scrna_integration")):
        raise RuntimeError(
            "找不到项目根目录（未检测到 src/scrna_integration/ 目录）。"
            "请从项目目录内启动 notebook。"
        )

    # ---- 2. 设置 BLAS 线程数 ----
    # A800 64核服务器上 OpenBLAS 默认全开 → 200+ 线程冻结。
    # 这里限制为 4，避免 scVI 等底层 BLAS 调用时线程爆炸。
    # 用 setdefault 而非直接赋值——不覆盖用户主动设置的值。
    _env_vars = {
        "OPENBLAS_NUM_THREADS": "4",
        "OMP_NUM_THREADS": "4",
        "MKL_NUM_THREADS": "4",
        "NUMBA_NUM_THREADS": "4",
    }
    for key, val in _env_vars.items():
        os.environ.setdefault(key, val)

    # ---- 3. 注入 sys.path ----
    _src_path = os.path.join(_root, "src")
    _already_in_path = _src_path in sys.path
    if not _already_in_path:
        sys.path.insert(0, _src_path)

    # ---- verbose：platform 决策可见 ----
    if verbose:
        print("=" * 50)
        print("Bootstrap 环境诊断")
        print("=" * 50)
        print(f"项目根路径 : {_root}")
        print(f"sys.path   : {'已存在' if _already_in_path else '已添加'} {_src_path}")
        print("BLAS 线程限制 :")
        for key, val in _env_vars.items():
            print(f"  {key} = {val}")
        print("=" * 50)

    return _root
