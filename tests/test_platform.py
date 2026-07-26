"""Tests for scrna_integration.platform — rscript_bin() path resolution."""

from __future__ import annotations

import os
import stat
import sys
import tempfile
from unittest.mock import MagicMock

import pytest

from scrna_integration.platform import (
    check_r_available, detect_device, env_check, platform_tag, rscript_bin,
)

# ---------------------------------------------------------------------------
# 辅助：在临时目录中模拟 conda 环境目录结构
# ---------------------------------------------------------------------------


def _make_fake_envs_dir(base: str, r_env_name: str = "scrna-integration-r") -> str:
    """在 base 下创建 envs/{r_env_name}/bin/Rscript 目录结构。

    返回 Rscript 的预期绝对路径。
    """
    envs_dir = os.path.join(base, "envs")
    r_env_dir = os.path.join(envs_dir, r_env_name)
    bin_dir = os.path.join(r_env_dir, "bin")
    os.makedirs(bin_dir, exist_ok=True)

    rscript_path = os.path.join(bin_dir, "Rscript")
    # 创建一个假的 Rscript 文件（内容无所谓，但要可执行）
    with open(rscript_path, "w") as f:
        f.write("#!/bin/bash\necho 'fake Rscript'\n")
    os.chmod(rscript_path, os.stat(rscript_path).st_mode | stat.S_IEXEC)

    return rscript_path


# ---------------------------------------------------------------------------
# 测试：从 CONDA_PREFIX 成功派生
# ---------------------------------------------------------------------------


def test_rscript_bin_from_conda_prefix(monkeypatch):
    """CONDA_PREFIX 已设置且派生路径存在 → 返回正确绝对路径。"""
    with tempfile.TemporaryDirectory() as tmpdir:
        expected = _make_fake_envs_dir(tmpdir)

        # 模拟当前激活的环境是 scrna-integration（Python 环境）
        current_env = os.path.join(tmpdir, "envs", "scrna-integration")
        os.makedirs(current_env, exist_ok=True)
        monkeypatch.setenv("CONDA_PREFIX", current_env)

        result = rscript_bin()
        assert result == expected
        assert os.path.isfile(result)


def test_rscript_bin_custom_env_name(monkeypatch):
    """传入自定义 r_env_name → 使用该环境名而非默认值。"""
    with tempfile.TemporaryDirectory() as tmpdir:
        expected = _make_fake_envs_dir(tmpdir, r_env_name="my-custom-r")

        current_env = os.path.join(tmpdir, "envs", "scrna-integration")
        os.makedirs(current_env, exist_ok=True)
        monkeypatch.setenv("CONDA_PREFIX", current_env)

        result = rscript_bin(r_env_name="my-custom-r")
        assert result == expected


# ---------------------------------------------------------------------------
# 测试：CONDA_PREFIX 已设置但派生路径不存在时回退到 shutil.which
# ---------------------------------------------------------------------------


def test_rscript_bin_fallback_to_which_when_derived_missing(monkeypatch, tmp_path):
    """CONDA_PREFIX 指向的位置下没有 scrna-integration-r → 回退到 PATH。"""
    # 防止 home 探测（优先级 1）意外命中本机真实 Rscript
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path)

    # 创建一个不存在的 envs 目录场景
    fake_prefix = str(tmp_path / "nonexistent_envs" / "scrna-integration")
    os.makedirs(fake_prefix, exist_ok=True)
    monkeypatch.setenv("CONDA_PREFIX", fake_prefix)

    # 在 tmp_path 下创建一个假的 Rscript，并把它加到 PATH
    fake_rscript = tmp_path / "Rscript"
    fake_rscript.write_text("#!/bin/bash\necho 'PATH Rscript'\n")
    fake_rscript.chmod(fake_rscript.stat().st_mode | stat.S_IEXEC)

    # 追加 PATH
    monkeypatch.setenv("PATH", str(tmp_path), prepend=os.pathsep)

    result = rscript_bin()
    assert result == str(fake_rscript)


# ---------------------------------------------------------------------------
# 测试：CONDA_PREFIX 未设置，但 home 目录下有常见 conda 安装
# ---------------------------------------------------------------------------


def test_rscript_bin_from_home_probing_when_no_conda_prefix(monkeypatch, tmp_path):
    """CONDA_PREFIX 未设置 + home 下有 miniforge3 → 命中探测。"""
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    # PATH 中不设 Rscript，确保只走 home 探测路径
    monkeypatch.setenv("PATH", str(tmp_path))

    # monkeypatch Path.home() 指向 tmp_path，在其中创建 miniforge3 结构
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path)
    home_install = tmp_path / "miniforge3"
    expected = _make_fake_envs_dir(str(home_install))

    result = rscript_bin()
    assert result == expected
    assert os.path.isfile(result)


# ---------------------------------------------------------------------------
# 测试：CONDA_PREFIX 未设置，但 PATH 中有 Rscript
# ---------------------------------------------------------------------------


def test_rscript_bin_from_path_when_no_conda_prefix(monkeypatch, tmp_path):
    """CONDA_PREFIX 未设置 → 回退到 shutil.which。"""
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    # 防止 home 探测（优先级 1）意外命中本机真实 Rscript
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path)

    fake_rscript = tmp_path / "Rscript"
    fake_rscript.write_text("#!/bin/bash\necho 'PATH only'\n")
    fake_rscript.chmod(fake_rscript.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", str(tmp_path), prepend=os.pathsep)

    result = rscript_bin()
    assert result == str(fake_rscript)


# ---------------------------------------------------------------------------
# 测试：所有方法都失败 → RuntimeError
# ---------------------------------------------------------------------------


def test_rscript_bin_raises_when_not_found(monkeypatch, tmp_path):
    """CONDA_PREFIX 未设置 + PATH 中无 Rscript → 抛出 RuntimeError。"""
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    # 防止 home 探测（优先级 1）意外命中本机真实 Rscript
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path)
    # 设置一个不包含 Rscript 的 PATH
    monkeypatch.setenv("PATH", str(tmp_path))

    with pytest.raises(RuntimeError, match="无法定位 Rscript"):
        rscript_bin()


def test_rscript_bin_raises_includes_env_name_in_message(monkeypatch, tmp_path):
    """异常消息中包含传入的 r_env_name。"""
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    # 防止 home 探测（优先级 1）意外命中本机真实 Rscript
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path)
    monkeypatch.setenv("PATH", str(tmp_path))

    with pytest.raises(RuntimeError, match="my-special-r"):
        rscript_bin(r_env_name="my-special-r")


# ---------------------------------------------------------------------------
# 测试：源代码不包含硬编码路径
# ---------------------------------------------------------------------------


def test_platform_module_no_hardcoded_paths():
    """platform.py 源码中不含 /Users/ 或 /home/ 硬编码绝对路径。"""
    import inspect

    import scrna_integration.platform as pm

    source = inspect.getsource(pm)
    # 允许在注释/文档字符串中出现这些路径作为示例说明
    # 但在实际代码行（非注释、非 docstring）中不应有
    lines = source.split("\n")
    in_docstring = False
    for lineno, line in enumerate(lines, start=1):
        stripped = line.strip()
        # 跳过空行
        if not stripped:
            continue
        # 跟踪多行文档字符串（三引号）的起止
        if '"""' in stripped or "'''" in stripped:
            # 计数三引号出现次数：奇数表示进入/退出 docstring
            triple_count = stripped.count('"""') + stripped.count("'''")
            if triple_count % 2 == 1:
                in_docstring = not in_docstring
            # 单行三引号（开头和结尾在同一行，如 """..."""）
            # triple_count 为偶数时不改变状态，但整行是 docstring
            if triple_count >= 2:
                continue
            if in_docstring or triple_count == 1:
                continue
        # 跳过注释行
        if stripped.startswith("#"):
            continue
        # 跳过 docstring 内部的中间行
        if in_docstring:
            continue
        # 检查是否包含硬编码绝对路径模式
        if "/Users/" in stripped or "/home/" in stripped:
            pytest.fail(
                f"platform.py 第 {lineno} 行包含硬编码路径: {stripped!r}"
            )


# ---------------------------------------------------------------------------
# 测试：shutil.which 返回的是绝对路径
# ---------------------------------------------------------------------------


def test_rscript_bin_returns_absolute_path(monkeypatch):
    """rscript_bin() 始终返回绝对路径。"""
    with tempfile.TemporaryDirectory() as tmpdir:
        expected = _make_fake_envs_dir(tmpdir)
        current_env = os.path.join(tmpdir, "envs", "scrna-integration")
        os.makedirs(current_env, exist_ok=True)
        monkeypatch.setenv("CONDA_PREFIX", current_env)

        result = rscript_bin()
        assert os.path.isabs(result)
        assert result == expected


# ---------------------------------------------------------------------------
# 测试：platform_tag() 平台标识映射
# ---------------------------------------------------------------------------


def test_platform_tag_linux_x86_64(monkeypatch):
    """Linux x86_64 → 'linux-64'。"""
    monkeypatch.setattr("platform.system", lambda: "Linux")
    monkeypatch.setattr("platform.machine", lambda: "x86_64")
    assert platform_tag() == "linux-64"


def test_platform_tag_darwin_arm64(monkeypatch):
    """macOS Apple Silicon → 'osx-arm64'。"""
    monkeypatch.setattr("platform.system", lambda: "Darwin")
    monkeypatch.setattr("platform.machine", lambda: "arm64")
    assert platform_tag() == "osx-arm64"


def test_platform_tag_darwin_x86_64(monkeypatch):
    """macOS Intel → 'osx-64'。"""
    monkeypatch.setattr("platform.system", lambda: "Darwin")
    monkeypatch.setattr("platform.machine", lambda: "x86_64")
    assert platform_tag() == "osx-64"


def test_platform_tag_unknown(monkeypatch):
    """无法识别的系统/架构 → '{system}-{machine}' 原样。"""
    monkeypatch.setattr("platform.system", lambda: "Windows")
    monkeypatch.setattr("platform.machine", lambda: "AMD64")
    assert platform_tag() == "windows-amd64"


def test_platform_tag_real_call():
    """真实调用不崩溃，返回值属已知平台标识之一（跨平台通用）。"""
    tag = platform_tag()
    assert isinstance(tag, str)
    assert len(tag) > 0
    # 跨平台断言：platform_tag() 在 Mac (osx-arm64/osx-64) 或 Linux (linux-64) 上
    # 均返回对应已知标识，不绑定单一平台
    assert tag in {"linux-64", "osx-arm64", "osx-64"}


# ---------------------------------------------------------------------------
# 测试：check_r_available() —— R 环境可用性检测（带标准化中文提示）
# ---------------------------------------------------------------------------


def test_check_r_available_success(monkeypatch):
    """mock rscript_bin 返回有效路径 → 返回 (path, True)。"""
    expected = "/path/to/Rscript"
    monkeypatch.setattr(
        "scrna_integration.platform.rscript_bin", lambda env="scrna-integration-r": expected
    )
    path, ok = check_r_available()
    assert ok is True
    assert path == expected


def test_check_r_available_failure(monkeypatch):
    """mock rscript_bin 抛出 RuntimeError → 返回 (None, False)。"""
    def _raise(*args, **kwargs):
        raise RuntimeError("无法定位 Rscript")
    monkeypatch.setattr("scrna_integration.platform.rscript_bin", _raise)
    path, ok = check_r_available()
    assert ok is False
    assert path is None


def test_check_r_available_prints_message(monkeypatch, capsys):
    """验证两种情况的打印输出均包含标准化中文提示。"""
    # 成功情况
    monkeypatch.setattr(
        "scrna_integration.platform.rscript_bin", lambda env="scrna-integration-r": "/fake/Rscript"
    )
    check_r_available()
    out = capsys.readouterr().out
    assert "R 环境就绪" in out

    # 失败情况
    def _raise(*args, **kwargs):
        raise RuntimeError("no")
    monkeypatch.setattr("scrna_integration.platform.rscript_bin", _raise)
    check_r_available()
    out = capsys.readouterr().out
    assert "R 环境未就绪" in out
    assert "conda env create -f environment-r.yml" in out


# ---------------------------------------------------------------------------
# 测试：env_check() —— 环境自检（返回结构 + verbose 不崩溃）
# ---------------------------------------------------------------------------


def test_env_check_returns_structure():
    """env_check 返回包含必需字段的 dict，不抛异常。"""
    result = env_check(verbose=False)
    assert isinstance(result, dict)
    for key in ["platform_tag", "conda_env", "ok", "checks", "warnings", "actions", "device"]:
        assert key in result
    assert isinstance(result["ok"], bool)
    assert isinstance(result["checks"], list)
    assert isinstance(result["actions"], list)


def test_env_check_no_crash_verbose():
    """verbose=True 不抛异常。"""
    env_check(verbose=True)  # 打印不应崩溃


# ---------------------------------------------------------------------------
# 测试：env_check() TF 冲突与 keras 残留检测
# ---------------------------------------------------------------------------


def test_env_check_tf_conflict_detected(monkeypatch):
    """主环境检测到 tensorflow 时应报 error（ok=False）。"""
    import importlib.metadata as _ilm
    from scrna_integration import platform as _plat
    _orig = _ilm.version

    def _fake_version(pkg):
        if pkg == "tensorflow":
            return "2.21.0"
        return _orig(pkg)

    monkeypatch.setattr(_ilm, "version", _fake_version)
    result = _plat.env_check(expected_env="scrna-integration", verbose=False)
    assert result["ok"] is False
    assert any("tensorflow" in w.lower() for w in result["warnings"])


def test_env_check_sccoda_env_no_tf_conflict(monkeypatch):
    """非主环境（如 scrna-sccoda）有 TF 不应报 error。"""
    import importlib.metadata as _ilm
    from scrna_integration import platform as _plat
    _orig = _ilm.version

    def _fake_version(pkg):
        if pkg == "tensorflow":
            return "2.21.0"
        return _orig(pkg)

    monkeypatch.setattr(_ilm, "version", _fake_version)
    result = _plat.env_check(expected_env="scrna-sccoda", verbose=False)
    # scCODA 环境本就该有 TF，不触发 TF 冲突 error
    assert not any("不应在主环境" in str(c) for c in result["checks"])


# ---------------------------------------------------------------------------
# 测试：detect_device() —— 计算设备自适应检测（ADR-0013）
# ---------------------------------------------------------------------------



def test_detect_device_cuda_available(monkeypatch):
    """cuda 可用 → accelerator=="gpu", device_str=="cuda"."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = True
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device()
    assert result["accelerator"] == "gpu"
    assert result["device_str"] == "cuda"
    assert result["devices"] == "auto"


def test_detect_device_mps_scvi_fallback_to_cpu(monkeypatch):
    """cuda=False/mps=True/for_method='scvi' → accelerator=='cpu'."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = False
    mock_torch.backends.mps.is_available.return_value = True
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device(prefer="auto", for_method="scvi")
    assert result["accelerator"] == "cpu"
    assert result["device_str"] == "cpu"
    assert "scVI" in result["reason"]


def test_detect_device_mps_no_method_uses_mps(monkeypatch):
    """cuda=False/mps=True/for_method=None → accelerator=='mps'."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = False
    mock_torch.backends.mps.is_available.return_value = True
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device(prefer="auto", for_method=None)
    assert result["accelerator"] == "mps"
    assert result["device_str"] == "mps"


def test_detect_device_neither_gpu_nor_mps(monkeypatch):
    """cuda=False/mps=False → accelerator=='cpu', device_str=='cpu'."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = False
    mock_torch.backends.mps.is_available.return_value = False
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device()
    assert result["accelerator"] == "cpu"
    assert result["device_str"] == "cpu"


def test_detect_device_torch_missing(monkeypatch):
    """torch 缺失（ImportError）→ cpu + reason 含 'torch 未安装'."""
    import sys
    from scrna_integration.platform import detect_device

    # 将 sys.modules 中的 torch 设为 None，模拟未安装
    monkeypatch.setitem(sys.modules, "torch", None)
    result = detect_device()
    assert result["accelerator"] == "cpu"
    assert result["device_str"] == "cpu"
    assert "torch 未安装" in result["reason"]


def test_detect_device_prefer_cpu_overrides_cuda(monkeypatch):
    """prefer='cpu' 即使 cuda 可用仍返回 cpu."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = True
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device(prefer="cpu")
    assert result["accelerator"] == "cpu"
    assert result["device_str"] == "cpu"
    assert "显式指定 CPU" in result["reason"]


def test_detect_device_prefer_cuda_unavailable_fallback(monkeypatch):
    """prefer='cuda' 但 cuda=False → 降级 cpu，reason 含警告."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = False
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device(prefer="cuda")
    assert result["accelerator"] == "cpu"
    assert result["device_str"] == "cpu"
    assert "不可用" in result["reason"]


def test_detect_device_prefer_cuda_available(monkeypatch):
    """prefer='cuda' + cuda 可用 → accelerator=='gpu'."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = True
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device(prefer="cuda")
    assert result["accelerator"] == "gpu"
    assert result["device_str"] == "cuda"
    assert "显式指定 CUDA" in result["reason"]


def test_detect_device_prefer_mps_available(monkeypatch):
    """prefer='mps' + mps 可用 → accelerator=='mps'."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = False
    mock_torch.backends.mps.is_available.return_value = True
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device(prefer="mps")
    assert result["accelerator"] == "mps"
    assert result["device_str"] == "mps"
    assert "显式指定 MPS" in result["reason"]


def test_detect_device_prefer_mps_unavailable_fallback(monkeypatch):
    """prefer='mps' 但 mps=False → 降级 cpu."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = False
    mock_torch.backends.mps.is_available.return_value = False
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device(prefer="mps")
    assert result["accelerator"] == "cpu"
    assert result["device_str"] == "cpu"
    assert "不可用" in result["reason"]


def test_detect_device_prefer_invalid_treated_as_auto(monkeypatch):
    """非法 prefer 值 → 按 auto 处理，reason 含注记."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = False
    mock_torch.backends.mps.is_available.return_value = False
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device(prefer="gpu")
    assert result["accelerator"] == "cpu"
    assert result["device_str"] == "cpu"
    assert "非法" in result["reason"]


def test_detect_device_mps_scanvi_fallback_to_cpu(monkeypatch):
    """cuda=False/mps=True/for_method='scanvi' → accelerator=='cpu'."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = False
    mock_torch.backends.mps.is_available.return_value = True
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device(prefer="auto", for_method="scanvi")
    assert result["accelerator"] == "cpu"
    assert result["device_str"] == "cpu"
    assert "scANVI" in result["reason"]


def test_detect_device_mps_sccraft_fallback_to_cpu(monkeypatch):
    """cuda=False/mps=True/for_method='sccraft' → accelerator=='cpu'."""
    from scrna_integration.platform import detect_device

    mock_torch = MagicMock()
    mock_torch.cuda.is_available.return_value = False
    mock_torch.backends.mps.is_available.return_value = True
    monkeypatch.setitem(sys.modules, "torch", mock_torch)

    result = detect_device(prefer="auto", for_method="sccraft")
    assert result["accelerator"] == "cpu"
    assert result["device_str"] == "cpu"
    assert "scCRAFT" in result["reason"]


def test_detect_device_return_keys():
    """返回 dict 含 accelerator/devices/device_str/reason 四字段."""
    from scrna_integration.platform import detect_device

    result = detect_device()
    for key in ["accelerator", "devices", "device_str", "reason"]:
        assert key in result, f"缺少字段: {key}"
    assert isinstance(result["accelerator"], str)
    assert isinstance(result["device_str"], str)
    assert isinstance(result["reason"], str)
    assert result["accelerator"] in ("gpu", "cpu", "mps")


# ---------------------------------------------------------------------------
# P1 补测：env_check Keras 残留 / 核心包缺失 / conda 不匹配
# ---------------------------------------------------------------------------


class TestEnvCheckExtended:
    """env_check 的补充检测分支——全用 mock 注入，不依赖真实环境。"""

    def test_env_check_keras_residual_detected(self, monkeypatch):
        """TF 已移除但 keras 残留 → warn。"""
        import importlib.metadata as _ilm
        from scrna_integration import platform as _plat
        _orig = _ilm.version

        def _fake_version(pkg):
            if pkg == "keras":
                return "3.9.0"
            if pkg == "tensorflow":
                return None  # TF 已移除
            return _orig(pkg)

        monkeypatch.setattr(_ilm, "version", _fake_version)
        result = _plat.env_check(expected_env="scrna-integration", verbose=False)
        # keras 残留检测记录在 checks 列表的 ("warn", msg) 元组中
        assert any(
            "keras" in str(c).lower()
            for c in result["checks"]
        )

    def test_env_check_core_pkg_missing_scanpy(self, monkeypatch):
        """核心包 scanpy 缺失 → ok=False + error 项。"""
        import importlib.metadata as _ilm
        from scrna_integration import platform as _plat
        _orig = _ilm.version

        def _fake_version(pkg):
            if pkg == "scanpy":
                return None  # 未安装
            return _orig(pkg)

        monkeypatch.setattr(_ilm, "version", _fake_version)
        result = _plat.env_check(expected_env="scrna-integration", verbose=False)
        assert result["ok"] is False
        assert any("scanpy" in str(c) for c in result["checks"])

    def test_env_check_conda_env_mismatch(self, monkeypatch):
        """当前 conda 环境名与期望不一致 → warn。"""
        from scrna_integration import platform as _plat
        monkeypatch.setenv("CONDA_DEFAULT_ENV", "base")
        monkeypatch.delenv("CONDA_PREFIX", raising=False)
        result = _plat.env_check(expected_env="scrna-integration", verbose=False)
        assert any("期望" in w for w in result["warnings"])


# ---------------------------------------------------------------------------
# P1 补测：detect_device for_method 未知值时行为
# ---------------------------------------------------------------------------


class TestDetectDeviceEdgeCases:
    """detect_device 边界场景。"""

    def test_detect_device_for_method_unknown_uses_mps(self, monkeypatch):
        """cuda=False/mps=True/for_method='unknown' → 以 MPS 处理（非 scvi/scanvi/sccraft）。"""
        from scrna_integration.platform import detect_device

        mock_torch = MagicMock()
        mock_torch.cuda.is_available.return_value = False
        mock_torch.backends.mps.is_available.return_value = True
        monkeypatch.setitem(sys.modules, "torch", mock_torch)

        result = detect_device(prefer="auto", for_method="unknown")
        # 非已知方法 → 不触发 scVI/scANVI/scCRAFT 的 CPU 回退 → 使用 MPS
        assert result["accelerator"] == "mps"
        assert result["device_str"] == "mps"
