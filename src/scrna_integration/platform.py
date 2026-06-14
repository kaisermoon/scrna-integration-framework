"""平台环境检测与路径解析（ADR-0010：OS 检测单点收口；ADR-0013：device 自适应层）。

本模块是框架中**唯一**允许感知操作系统/平台差异的位置。
所有其他模块（src/ 其余代码、notebooks/）**禁止**直接使用
``sys.platform``、``os.uname`` 或硬编码平台绝对路径（``/Users/``、
``/home/`` 等）；需要平台相关路径时统一通过本模块获取。

自 ADR-0013 起，计算设备检测（CUDA / MPS / CPU）也收口至本模块的
:func:`detect_device`，notebook/src 其他位置不再自行判断
``torch.cuda.is_available()`` 或 ``torch.backends.mps``。

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


def check_r_available(r_env_name: str = "scrna-integration-r") -> tuple:
    """检查 R 环境是否可用。

    内部调用 :func:`rscript_bin`，捕获 ``RuntimeError``。
    成功返回 ``(path, True)``，失败返回 ``(None, False)`` 并打印标准化提示。

    Parameters
    ----------
    r_env_name : str
        conda R 环境名称，默认 ``"scrna-integration-r"``。

    Returns
    -------
    tuple[str | None, bool]
        ``(Rscript 绝对路径, True)`` 或 ``(None, False)``。

    Examples
    --------
    >>> from scrna_integration.platform import check_r_available
    >>> RSCRIPT_BIN, R_AVAILABLE = check_r_available()
    >>> if R_AVAILABLE:
    ...     # 调用 R 工具
    ...     subprocess.run([RSCRIPT_BIN, "--vanilla", script])
    """
    try:
        path = rscript_bin(r_env_name)
    except RuntimeError:
        print(
            f"⚠️ R 环境未就绪（{r_env_name}），"
            f"R 依赖分析将跳过。如需启用，请安装 conda 环境: "
            f"conda env create -f environment-r.yml"
        )
        return (None, False)
    print(f"✓ R 环境就绪: {path}")
    return (path, True)


def detect_device(prefer="auto", for_method=None):
    """检测当前计算设备并返回训练配置参数。

    根据 ADR-0013，统一收口设备检测逻辑，支持 Mac-MPS / Linux-CPU /
    Ubuntu-CUDA 三环境自适应。返回的 dict 可直接解包喂给
    ``scvi.model.SCVI.train(accelerator=..., devices=...)``。

    检测优先级：
    1. try import torch（局部 import）；ImportError → 回退 CPU
    2. CUDA 可用 > MPS 可用 > CPU
    3. 显式 ``prefer`` 参数可覆盖自动检测
    4. ``for_method`` 允许方法级默认行为（scVI/scANVI 在 Mac 上默认
       走 CPU 而非 MPS，因 MPS 数值稳定性未经验证）

    参数
    ----------
    prefer : str
        显式设备偏好。合法值：
        - ``"auto"``（默认）——自动检测
        - ``"cuda"``——强制 CUDA GPU
        - ``"mps"``——强制 Apple MPS
        - ``"cpu"``——强制 CPU
        非法值视为 ``"auto"`` 并在 reason 中注明。
    for_method : str or None
        调用方方法名（``"scvi"`` / ``"scanvi"`` / ``"sccraft"``），
        用于 auto 模式下的方法级默认行为。None 表示不区分方法。

    返回
    -------
    dict
        - ``"accelerator"`` (str) —— ``"gpu"`` / ``"cpu"`` / ``"mps"``，
          直接传给 ``model.train(accelerator=...)``
        - ``"devices"`` (object) —— ``"auto"`` / ``1`` / ``[0]``，
          直接传给 ``model.train(devices=...)``
        - ``"device_str"`` (str) —— ``"cuda"`` / ``"mps"`` / ``"cpu"``
          （torch 风格，日志用）
        - ``"reason"`` (str) —— 中文决策说明

    使用示例
    --------
    >>> from scrna_integration.platform import detect_device
    >>> _dev = detect_device("auto", for_method="scvi")
    >>> model.train(accelerator=_dev["accelerator"], devices=_dev["devices"])
    """
    # 1) 局部 import torch（支持测试 monkeypatch）
    try:
        import torch
    except ImportError:
        return {
            "accelerator": "cpu",
            "devices": "auto",
            "device_str": "cpu",
            "reason": "torch 未安装，回退 CPU",
        }

    cuda_ok = torch.cuda.is_available()
    # 注意：torch.backends.mps 仅在 macOS PyTorch >= 1.12 存在
    mps_ok = hasattr(torch.backends, "mps") and torch.backends.mps.is_available()

    # --- 合法值校验 ---
    _valid = {"auto", "cuda", "mps", "cpu"}
    _prefer = prefer
    if _prefer not in _valid:
        _prefer = "auto"

    # 2) 显式 prefer 覆盖
    if _prefer == "cuda":
        if cuda_ok:
            return {"accelerator": "gpu", "devices": "auto",
                    "device_str": "cuda", "reason": "显式指定 CUDA"}
        else:
            return {"accelerator": "cpu", "devices": "auto",
                    "device_str": "cpu",
                    "reason": "请求 CUDA 但不可用，回退 CPU"}

    if _prefer == "mps":
        if mps_ok:
            return {"accelerator": "mps", "devices": "auto",
                    "device_str": "mps", "reason": "显式指定 MPS"}
        else:
            return {"accelerator": "cpu", "devices": "auto",
                    "device_str": "cpu",
                    "reason": "请求 MPS 但不可用，回退 CPU"}

    if _prefer == "cpu":
        return {"accelerator": "cpu", "devices": "auto",
                "device_str": "cpu",
                "reason": "显式指定 CPU"}

    # 3) auto 模式
    # 非法 prefer 已在上面被转为 auto，此处只追加注记
    _note = f"（prefer={prefer!r} 非法，已按 auto 处理）" if prefer not in _valid else ""

    if cuda_ok:
        return {"accelerator": "gpu", "devices": "auto",
                "device_str": "cuda",
                "reason": "检测到 CUDA GPU" + _note}

    if mps_ok:
        if for_method in ("scvi", "scanvi"):
            return {
                "accelerator": "cpu", "devices": "auto",
                "device_str": "cpu",
                "reason": (
                    "Mac/MPS：scVI/scANVI 默认 CPU"
                    "（MPS 数值稳定性未验证），如需 MPS 请显式 DEVICE='mps'"
                ) + _note,
            }
        if for_method == "sccraft":
            return {
                "accelerator": "cpu", "devices": "auto",
                "device_str": "cpu",
                "reason": "scCRAFT 内部硬编码 CPU，device 不可控（见 ADR-0013）" + _note,
            }
        # for_method is None → MPS OK
        return {"accelerator": "mps", "devices": "auto",
                "device_str": "mps",
                "reason": "Mac 检测到 MPS（具体方法可能回退 CPU）" + _note}

    # 4) 无 GPU / 无 MPS
    return {"accelerator": "cpu", "devices": "auto",
            "device_str": "cpu",
            "reason": "未检测到 GPU/MPS，使用 CPU" + _note}


def env_check(expected_env="scrna-integration", verbose=True, device_prefer="auto"):
    """诊断当前运行环境是否符合本项目期望。

    检测三层：(1) 平台标识 (2) 激活的 conda 环境 (3) 关键包存在性与已知冲突
    (4) 计算设备自适应检测。遵循 ADR-0010 哲学：只诊断不修改环境——对齐决策归人。
    给出需要的调整命令，但绝不自动执行。

    设计为可在每个 notebook 的 setup cell 一行调用：
        from scrna_integration.platform import env_check
        env_check()

    Parameters
    ----------
    expected_env : str
        本 stage 期望的 conda 环境名。主流水线用 "scrna-integration"。
    verbose : bool
        是否打印诊断报告（默认 True）。
    device_prefer : str
        传给 :func:`detect_device` 的设备偏好（默认 ``"auto"``）。
        当 notebook 显式指定了 ``DEVICE="cpu"`` 等值时，传入相同值可使
        env_check 显示的设备与后续训练保持一致。

    Returns
    -------
    dict
        {platform_tag, conda_env, ok(bool), checks(list), warnings(list),
         actions(list), device(str|None)}
    """
    import importlib.metadata as _ilm
    import importlib.util as _ilu

    _tag = platform_tag()
    # 当前激活的 conda 环境名
    _conda_env = os.environ.get("CONDA_DEFAULT_ENV", "")
    if not _conda_env:
        _prefix = os.environ.get("CONDA_PREFIX", "")
        _conda_env = os.path.basename(_prefix) if _prefix else "unknown"

    _checks = []    # [(level, msg)]  level in {"ok","warn","error"}
    _warnings = []
    _actions = []   # 需要 PI 执行的调整命令

    # --- 平台 ---
    _checks.append(("ok", f"平台: {_tag}"))

    # --- conda 环境匹配 ---
    if _conda_env == expected_env:
        _checks.append(("ok", f"conda 环境: {_conda_env}"))
    else:
        _checks.append(("warn", f"conda 环境: {_conda_env}（期望 {expected_env}）"))
        _warnings.append(f"当前 conda 环境是 '{_conda_env}'，本 stage 期望 '{expected_env}'")
        _actions.append(f"conda activate {expected_env}")

    # --- 关键包检测 ---
    def _ver(pkg):
        try:
            return _ilm.version(pkg)
        except _ilm.PackageNotFoundError:
            return None

    # 必需的核心包
    for _pkg in ["scanpy", "anndata", "numpy"]:
        _v = _ver(_pkg)
        if _v:
            _checks.append(("ok", f"{_pkg}={_v}"))
        else:
            _checks.append(("error", f"{_pkg} 未安装（核心依赖缺失）"))
            _warnings.append(f"核心包 {_pkg} 缺失——环境可能不对")

    # 嵌入/注释 stage 需要的包（缺失只警告不报错——不是所有 stage 都用）
    for _pkg in ["torch", "scvi-tools"]:
        _v = _ver(_pkg)
        if _v:
            _checks.append(("ok", f"{_pkg}={_v}"))
        else:
            _checks.append(("warn", f"{_pkg} 未安装（scVI/scANVI/scCRAFT 步骤需要）"))

    # --- 计算设备检测（复用 detect_device，仅诊断）---
    _device_str = None
    try:
        _dev = detect_device(prefer=device_prefer)
        _checks.append(("ok", f"计算设备: {_dev['device_str']}（{_dev['reason']}）"))
        _device_str = _dev["device_str"]
    except Exception as _e:
        _checks.append(("warn", f"设备检测跳过: {_e}"))

    # --- TF 冲突检测（主环境的关键检查）---
    if expected_env == "scrna-integration":
        _tf = _ver("tensorflow")
        _keras = _ver("keras")
        if _tf:
            _checks.append(("error", f"tensorflow={_tf} 不应在主环境！"))
            _warnings.append(
                "检测到 tensorflow 在主环境——它的 oneDNN/MKL 与 PyTorch backward "
                "在同进程冲突会导致 scCRAFT/scVI 训练 segfault（见 ADR-0012）"
            )
            _actions.append("pip uninstall -y tensorflow tensorflow-probability keras  "
                            "# scCODA 请用 scrna-sccoda 环境")
        elif _keras:
            # keras 3.x 可独立于 TF，但残留可能引入混淆/依赖——提示清理
            _checks.append(("warn", f"keras={_keras} 残留（tensorflow 已移除，建议一并清理）"))
            _actions.append("pip uninstall -y keras")

    # --- scCRAFT 可选提示 ---
    _sccraft_installed = False
    try:
        _sccraft_installed = _ilu.find_spec("scCRAFT") is not None
    except Exception:
        pass
    if _sccraft_installed:
        _checks.append(("ok", "scCRAFT 已安装（04 可用 EMBEDDING_METHODS 加 'sccraft'）"))
    else:
        _checks.append(("warn", "scCRAFT 未安装（如需用 scCRAFT 嵌入: "
                        "git clone https://github.com/ch2343/scCRAFT && pip install . "
                        "&& pip uninstall -y tensorflow tensorflow-probability keras）"))

    # --- 平台同步提示（Mac 首次运行）---
    if _tag.startswith("osx"):
        _checks.append(("warn", "macOS 平台——如首次在 Mac 运行或环境久未同步，请查看 docs/MAC-SYNC.md"))
        _warnings.append("Mac 环境同步清单见 docs/MAC-SYNC.md（本项目最近在 Linux 做过环境调整）")

    _ok = not any(level == "error" for level, _ in _checks)

    if verbose:
        print("=" * 56)
        print(f"环境检测  |  平台 {_tag}  |  conda env: {_conda_env}")
        print("=" * 56)
        _icon = {"ok": "✓", "warn": "⚠️", "error": "❌"}
        for _level, _msg in _checks:
            print(f"  {_icon[_level]} {_msg}")
        if _actions:
            print("\n  需要的调整（请手动执行，本函数不自动改环境）:")
            for _a in _actions:
                print(f"    $ {_a}")
        if _ok:
            if _warnings:
                print("\n  状态: ⚠️ 可运行但有提示（见上）")
            else:
                print("\n  状态: ✓ 环境就绪")
        else:
            print("\n  状态: ❌ 环境有问题，请先处理上方 error 项")
        print("=" * 56)

    return {
        "platform_tag": _tag,
        "conda_env": _conda_env,
        "ok": _ok,
        "checks": _checks,
        "warnings": _warnings,
        "actions": _actions,
        "device": _device_str,
    }
