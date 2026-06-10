"""Tests for scrna_integration.platform — rscript_bin() path resolution."""

from __future__ import annotations

import os
import stat
import tempfile

import pytest

from scrna_integration.platform import platform_tag, rscript_bin

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
# 测试：CONDA_PREFIX 未设置，但 PATH 中有 Rscript
# ---------------------------------------------------------------------------


def test_rscript_bin_from_path_when_no_conda_prefix(monkeypatch, tmp_path):
    """CONDA_PREFIX 未设置 → 回退到 shutil.which。"""
    monkeypatch.delenv("CONDA_PREFIX", raising=False)

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
    # 设置一个不包含 Rscript 的 PATH
    monkeypatch.setenv("PATH", str(tmp_path))

    with pytest.raises(RuntimeError, match="无法定位 Rscript"):
        rscript_bin()


def test_rscript_bin_raises_includes_env_name_in_message(monkeypatch, tmp_path):
    """异常消息中包含传入的 r_env_name。"""
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
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
    """真实调用不崩溃，返回非空字符串。"""
    tag = platform_tag()
    assert isinstance(tag, str)
    assert len(tag) > 0
    # 本机是 Alibaba Cloud Linux 3 x86_64，应返回 "linux-64"
    assert tag == "linux-64"
