"""平台环境检测与路径解析（ADR-0010：OS 检测单点收口）。

本模块是框架中**唯一**允许感知操作系统/平台差异的位置。
所有其他模块（src/ 其余代码、notebooks/）**禁止**直接使用
``sys.platform``、``os.uname`` 或硬编码平台绝对路径（``/Users/``、
``/home/`` 等）；需要平台相关路径时统一通过本模块获取。

为什么需要这个模块？
--------------------
框架的 conda 环境结构在 Mac（osx-arm64）和 Linux 服务器（linux-64）
上保持一致（通过精确 pin 的源 spec + env_parity 诊断脚本人工对齐保证），但 conda 安装根目录
在不同机器上不同。例如：
  - Mac:  /Users/alice/miniforge3/envs/scrna-integration-r/bin/Rscript
  - Linux: /home/bob/miniforge3/envs/scrna-integration-r/bin/Rscript
  - 本机（巧合）：/Users/zhongzishao/miniforge3/...（Alibaba Cloud Linux，
    用户名恰好是 zhongzishao，路径含 /Users/ 纯属巧合，不是 Mac）

如果 notebook 里硬编码绝对路径，换一台机器就坏。正确做法是从运行时
环境变量动态推导——conda 激活后会设置 ``CONDA_PREFIX`` 指向当前环境
的安装目录，同 conda 安装下的其他环境路径可以从它派生出来。

回退链设计
----------
1. 首选：从 CONDA_PREFIX 环境变量派生（conda 环境激活时自动设置）。
   推导公式：``dirname(CONDA_PREFIX) / {r_env_name} / bin / Rscript``。
   这覆盖了 99% 的正常使用场景。
2. 回退一：调用 ``shutil.which("Rscript")`` 在系统 PATH 中查找。
   用于 conda 未激活但 Rscript 已在 PATH 中的场景。
3. 回退二（最终兜底）：抛出清晰的中文异常，提示用户激活或创建
   ``scrna-integration-r`` 环境。**绝不**写死任何机器特定的绝对路径，
   因为"写死在我机器上能跑"就是"换台机器必坏"。
"""

from __future__ import annotations

import os
import platform as _platform
import shutil


def platform_tag() -> str:
    """返回标准化的平台标识字符串，用于快照文件命名与跨机器对比。

    为什么需要这个函数？
    --------------------
    项目同时跑在 Mac（Apple Silicon）和 Linux 服务器（x86_64）上。
    两台机器的 conda 环境包清单需要分别导出、统一对比。
    如果没有一个标准化的平台标识，对比脚本就不知道"谁和谁比"。
    本函数是 ADR-0010「OS 检测单点收口」的合理扩充——平台标识也属 OS 感知。

    返回
    -------
    str
        标准化平台标识：
        - ``"linux-64"``   —— Linux + x86_64
        - ``"osx-arm64"``  —— macOS + Apple Silicon（arm64）
        - ``"osx-64"``     —— macOS + Intel x86_64
        - 其他无法识别时返回 ``"{system}-{machine}"`` 原样（如 ``"windows-AMD64"``）

    使用示例
    --------
    >>> from scrna_integration.platform import platform_tag
    >>> tag = platform_tag()   # 本机返回 "linux-64"
    >>> snapshot_path = f"docs/env-snapshots/{tag}.json"
    """
    system = _platform.system().lower()    # "linux" / "darwin" / "windows"
    machine = _platform.machine().lower()  # "x86_64" / "arm64" / "amd64"

    # macOS 需要区分 Intel 和 Apple Silicon
    if system == "darwin":
        if machine == "arm64":
            return "osx-arm64"
        if machine in ("x86_64", "amd64"):
            return "osx-64"
        return f"osx-{machine}"

    # Linux：绝大部分情况是 x86_64 服务器
    if system == "linux":
        if machine == "x86_64":
            return "linux-64"
        # 将来的 aarch64 服务器场景
        if machine == "aarch64":
            return "linux-aarch64"
        return f"linux-{machine}"

    # 无法识别时返回原始 system-machine 组合，不造缩写
    return f"{system}-{machine}"


def rscript_bin(r_env_name: str = "scrna-integration-r") -> str:
    """返回 conda R 环境中 Rscript 可执行文件的绝对路径。

    参数
    ----------
    r_env_name : str
        conda R 环境名称，默认 ``"scrna-integration-r"``。
        如果你的 R 环境名不同（比如 ``"my-r-env"``），传入即可。

    返回
    -------
    str
        Rscript 可执行文件的绝对路径，例如：
        ``/home/bob/miniforge3/envs/scrna-integration-r/bin/Rscript``

    异常
    -------
    RuntimeError
        当 CONDA_PREFIX 未设置、派生路径不存在、且 PATH 中也找不到
        Rscript 时抛出。异常消息包含明确的修复指引。

    使用示例
    --------
    >>> from scrna_integration.platform import rscript_bin
    >>> RSCRIPT_BIN = rscript_bin()
    >>> # 如果环境名不同：
    >>> RSCRIPT_BIN = rscript_bin(r_env_name="scrna-integration-r-dev")
    """
    # ---- 首选：从当前激活的 conda 环境推导 ----
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        # conda 环境目录结构：
        #   {install_root}/envs/{env_name}/bin/Rscript
        # CONDA_PREFIX 指向当前环境：{install_root}/envs/{current_env}
        # 取上一级得到 {install_root}/envs/，再拼 R 环境名
        envs_dir = os.path.dirname(conda_prefix)
        candidate = os.path.join(envs_dir, r_env_name, "bin", "Rscript")
        if os.path.isfile(candidate):
            return candidate

    # ---- 回退一：在系统 PATH 中查找 Rscript ----
    which_result = shutil.which("Rscript")
    if which_result is not None:
        return which_result

    # ---- 回退二：无法定位，抛出清晰异常 ----
    raise RuntimeError(
        f"无法定位 Rscript 可执行文件。\n"
        f"\n"
        f"已尝试以下方式：\n"
        f"  1. 从 CONDA_PREFIX 派生同 conda 安装下的 "
        f"'{r_env_name}' 环境路径\n"
        f"     → 当前 CONDA_PREFIX="
        f"{conda_prefix!r}（{'未设置' if not conda_prefix else '已设置，但派生路径不存在'}）\n"
        f"  2. 在系统 PATH 中搜索 'Rscript'\n"
        f"     → 未找到\n"
        f"\n"
        f"请选择以下任一方式修复：\n"
        f"  a) 激活 conda R 环境后重跑 notebook：\n"
        f"     conda activate {r_env_name}\n"
        f"  b) 如果环境尚未创建，使用项目提供的源 spec 创建：\n"
        f"     conda env create -f environment-r.yml\n"
        f"  c) 如果环境名不是默认的 '{r_env_name}'，请传入正确的名称：\n"
        f"     from scrna_integration.platform import rscript_bin\n"
        f"     RSCRIPT_BIN = rscript_bin(r_env_name='你的环境名')\n"
        f"  d) 如果 Rscript 已在 PATH 中但上述仍失败，请检查 "
        f"shutil.which('Rscript') 的返回值并报告。"
    )
