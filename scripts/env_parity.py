#!/usr/bin/env python
"""跨机器环境对齐诊断脚本（ADR-0010 修订版）。

为什么需要这个脚本？
--------------------
项目同时跑在 Mac（osx-arm64）和 Linux 服务器（linux-64）上。
两台机器各自按 environment.yml / environment-r.yml 装包（精确 pin），
但 conda 求解过程仍可能产出微小差异（传递依赖版本、channel 选择等）。
本脚本提供两种能力：
  1. snapshot —— 抓取当前机器 conda 环境包清单，导出为可对比的快照
  2. compare  —— 对比两份快照，列出差异供人工决策对齐

脚本不自动修改任何环境（硬约束：对齐决策归人，不外包给脚本）。

用法：
  python scripts/env_parity.py snapshot
  python scripts/env_parity.py compare [--a FILE] [--b FILE]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import socket
import subprocess
import sys
import time
from typing import Any

# ============================================================================
# 常量
# ============================================================================
# 项目根目录（scripts/ 的上一级）
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# ---- 自举导入 platform_tag（importlib 直加载源文件，不触发包 __init__） ----
# 不能 ``from scrna_integration.platform import platform_tag``，
# 因为 scrna_integration/__init__.py 触发的级联导入（anndata / scanpy 等）
# 只在 conda 环境内可用。改用 importlib.util 直接加载 platform.py 源文件，
# 不经过包 __init__，保持 ADR-0010 单点收口。
_spec = importlib.util.spec_from_file_location(
    "scrna_integration_platform",
    os.path.join(PROJECT_ROOT, "src", "scrna_integration", "platform.py"),
)
_platform_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_platform_mod)
platform_tag = _platform_mod.platform_tag

# 快照输出目录
SNAPSHOT_DIR = os.path.join(PROJECT_ROOT, "docs", "env-snapshots")

# 需要导出包清单的两个 conda 环境名
PY_ENV = "scrna-integration"
R_ENV = "scrna-integration-r"

# 默认对比的两台机器快照文件（相对于项目根）
DEFAULT_A = os.path.join(SNAPSHOT_DIR, "linux-64.json")
DEFAULT_B = os.path.join(SNAPSHOT_DIR, "osx-arm64.json")

# ============================================================================
# 工具函数：conda 命令封装
# ============================================================================


def run_conda(args: list[str], timeout: int = 60) -> subprocess.CompletedProcess:
    """安全地跑 conda 命令（stdout 捕获为文本，stderr 管道丢弃）。

    为什么不用 shell=True？
    -----------------------
    shell=True 带来注入风险和 shell 转义复杂性。直接传 list 给 subprocess.run
    是最干净的调用方式。conda 命令不需要 shell 解释（参数不含通配符/管道）。
    """
    return subprocess.run(
        ["conda"] + args,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def conda_list_json(env_name: str) -> list[dict[str, str]] | None:
    """对指定 conda 环境跑 ``conda list -n {env} --json`` 并返回包列表。

    为什么解析 JSON 而不是解析纯文本？
    ---------------------------------
    ``conda list --json`` 是 conda 官方提供的机器可读输出，不会因列对齐、
    空行、表头等文本格式变化而解析失败。JSON 天然适合脚本消费。

    返回
    -------
    list[dict] 或 None
        每个 dict 含 ``name`` / ``version`` / ``channel`` 等键。
        如果环境不存在（conda 命令失败）则返回 None。
    """
    result = run_conda(["list", "-n", env_name, "--json"], timeout=30)
    if result.returncode != 0:
        # conda list 失败通常意味着环境不存在
        print(f"  [警告] 环境 '{env_name}' 不存在或无法读取（conda list 返回非零退出码）")
        return None
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError:
        print(f"  [错误] 环境 '{env_name}' 的输出不是有效 JSON，请检查 conda 版本")
        return None


def build_package_dict(pkg_list: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """将 conda list JSON 输出转成 {包名: {version, channel}} 格式。

    为什么转成 dict 而不是保留 list？
    ---------------------------------
    快照对比需要按包名快速查找（O(1) 而非 O(n)），且同名包在不同环境下
    版本/来源可能不同。保留 version 和 channel 两个字段以便诊断差异根因。
    """
    result: dict[str, dict[str, str]] = {}
    for pkg in pkg_list:
        name = pkg.get("name", "")
        if not name:
            continue
        result[name] = {
            "version": pkg.get("version", "?"),
            "channel": pkg.get("channel", "?"),
        }
    return result


def check_gpu_available() -> bool:
    """检测本机是否有可用的 NVIDIA GPU。

    通过检查 ``nvidia-smi`` 命令是否可执行来判断。
    不实际调用 nvidia-smi（避免引入不必要的子进程延迟），
    只检查命令是否在 PATH 中。
    """
    return shutil.which("nvidia-smi") is not None


def conda_version() -> str:
    """获取 conda 自身版本号。"""
    result = run_conda(["--version"], timeout=10)
    if result.returncode == 0:
        return result.stdout.strip()
    return "unknown"


# ============================================================================
# snapshot 子命令：导出当前机器环境快照
# ============================================================================


def cmd_snapshot() -> int:
    """感知当前机器、导出 conda 环境包清单、写入 JSON 快照。

    退出码：0（成功）/ 1（部分环境读取失败，但快照仍写出）。
    """
    # ---- 1. 收集机器身份信息 ----
    tag = platform_tag()
    hostname = socket.gethostname()
    system_full = f"{os.uname().sysname} {os.uname().release}"
    python_version = sys.version
    cv = conda_version()
    has_gpu = check_gpu_available()

    print("=== 环境快照：当前机器身份 ===")
    print(f"  platform_tag:  {tag}")
    print(f"  hostname:      {hostname}")
    print(f"  system:        {system_full}")
    print(f"  python:        {python_version.split()[0]}")
    print(f"  conda:         {cv}")
    print(f"  GPU:           {'是' if has_gpu else '无（nvidia-smi 未检测到）'}")
    print()

    # ---- 2. 导出两个环境的包清单 ----
    had_errors = False

    print(f"--- 环境: {PY_ENV} ---")
    py_list = conda_list_json(PY_ENV)
    py_packages = build_package_dict(py_list) if py_list is not None else None
    if py_packages is not None:
        print(f"  包数: {len(py_packages)}")
    else:
        had_errors = True

    print(f"--- 环境: {R_ENV} ---")
    r_list = conda_list_json(R_ENV)
    r_packages = build_package_dict(r_list) if r_list is not None else None
    if r_packages is not None:
        print(f"  包数: {len(r_packages)}")
    else:
        had_errors = True

    # ---- 3. 组装快照 JSON ----
    snapshot: dict[str, Any] = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "generated_by": f"env_parity.py snapshot on {hostname}",
        "machine": {
            "platform_tag": tag,
            "hostname": hostname,
            "system": system_full,
            "python_version": python_version,
            "conda_version": cv,
            "gpu_available": has_gpu,
        },
        "environments": {
            PY_ENV: py_packages,
            R_ENV: r_packages,
        },
    }

    # ---- 4. 写入快照文件 ----
    os.makedirs(SNAPSHOT_DIR, exist_ok=True)
    out_path = os.path.join(SNAPSHOT_DIR, f"{tag}.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(snapshot, f, indent=2, ensure_ascii=False)
        f.write("\n")  # 末尾换行，POSIX 友好

    print(f"\n快照已写入: {out_path}")
    print("提交到 git 后，另一台机器跑 snapshot 即可用 compare 对比差异。")

    return 1 if had_errors else 0


# ============================================================================
# compare 子命令：对比两份快照
# ============================================================================


def _compare_one_env(
    env_name: str,
    pkgs_a: dict[str, dict[str, str]] | None,
    pkgs_b: dict[str, dict[str, str]] | None,
    label_a: str,
    label_b: str,
) -> int:
    """对比单个环境在两份快照中的差异，打印表格。

    返回该环境的差异计数（版本不一致 + 缺失包总数）。
    """
    if pkgs_a is None and pkgs_b is None:
        print(f"  [跳过] 两边环境 '{env_name}' 都不存在，无法对比")
        return 0
    if pkgs_a is None:
        print(f"  [跳过] 环境 '{env_name}' 仅在 {label_b} 存在（{len(pkgs_b)} 包），"
              f"{label_a} 快照中缺失")
        return len(pkgs_b) if pkgs_b else 0
    if pkgs_b is None:
        print(f"  [跳过] 环境 '{env_name}' 仅在 {label_a} 存在（{len(pkgs_a)} 包），"
              f"{label_b} 快照中缺失")
        return len(pkgs_a) if pkgs_a else 0

    # 收集差异
    all_names = set(pkgs_a.keys()) | set(pkgs_b.keys())
    version_diffs: list[tuple[str, str, str]] = []  # (pkg, ver_a, ver_b)
    only_a: list[str] = []
    only_b: list[str] = []

    for name in sorted(all_names):
        in_a = name in pkgs_a
        in_b = name in pkgs_b
        if in_a and in_b:
            if pkgs_a[name]["version"] != pkgs_b[name]["version"]:
                version_diffs.append(
                    (name, pkgs_a[name]["version"], pkgs_b[name]["version"])
                )
        elif in_a and not in_b:
            only_a.append(name)
        else:
            only_b.append(name)

    diff_count = len(version_diffs) + len(only_a) + len(only_b)

    # 打印该环境标题
    print(f"\n{'=' * 70}")
    print(f"  环境: {env_name}")
    print(f"  文件 A ({label_a}): {len(pkgs_a)} 包")
    print(f"  文件 B ({label_b}): {len(pkgs_b)} 包")
    print(f"  差异项: {diff_count}")
    print(f"{'=' * 70}")

    if diff_count == 0:
        print("  ✓ 完全一致")
        return 0

    # 1. 版本不一致
    if version_diffs:
        print(f"\n  【版本不一致】（{len(version_diffs)} 项）")
        # 列宽自适应：包名最长的那一行
        max_name = max(len(row[0]) for row in version_diffs)
        max_ver_a = max(len(row[1]) for row in version_diffs)
        max_ver_b = max(len(row[2]) for row in version_diffs)
        header = (
            f"  {'包名':<{max_name}}  "
            f"{'A版本 (' + label_a + ')':<{max_ver_a + 4}}  "
            f"{'B版本 (' + label_b + ')':<{max_ver_b + 4}}"
        )
        print(header)
        print(f"  {'-' * (max_name + max_ver_a + max_ver_b + 15)}")
        for name, va, vb in version_diffs:
            print(f"  {name:<{max_name}}  {va:<{max_ver_a + 4}}  {vb:<{max_ver_b + 4}}")

    # 2. 仅 A 有
    if only_a:
        print(f"\n  【仅 {label_a} 有】（{len(only_a)} 项）")
        for name in only_a:
            va = pkgs_a[name]["version"]
            ch = pkgs_a[name]["channel"]
            print(f"    {name}=={va}  (channel: {ch})")

    # 3. 仅 B 有
    if only_b:
        print(f"\n  【仅 {label_b} 有】（{len(only_b)} 项）")
        for name in only_b:
            vb = pkgs_b[name]["version"]
            ch = pkgs_b[name]["channel"]
            print(f"    {name}=={vb}  (channel: {ch})")

    return diff_count


def cmd_compare(args: argparse.Namespace) -> int:
    """对比两份环境快照，输出差异清单。

    退出码：0（完全一致）/ 1（有差异）/ 2（快照文件缺失，无法对比）。
    """
    path_a = args.a or DEFAULT_A
    path_b = args.b or DEFAULT_B

    # 检查快照文件是否存在
    missing = []
    for label, path in [("A", path_a), ("B", path_b)]:
        if not os.path.isfile(path):
            missing.append((label, path))

    if missing:
        for label, path in missing:
            basename = os.path.basename(path)
            print(f"[缺失] 快照文件 {label} 不存在: {path}")
            print("        请在对应机器上运行：python scripts/env_parity.py snapshot")
            print(f"        生成 {basename} 后提交到 git，再跑 compare。")
        return 2

    # 加载快照
    try:
        with open(path_a, encoding="utf-8") as f:
            snap_a = json.load(f)
    except json.JSONDecodeError as e:
        print(f"[错误] 快照文件 A 不是有效 JSON: {path_a}")
        print(f"       {e}")
        return 2
    try:
        with open(path_b, encoding="utf-8") as f:
            snap_b = json.load(f)
    except json.JSONDecodeError as e:
        print(f"[错误] 快照文件 B 不是有效 JSON: {path_b}")
        print(f"       {e}")
        return 2

    # 提取身份信息用于显示
    machine_a = snap_a.get("machine", {})
    machine_b = snap_b.get("machine", {})
    label_a = machine_a.get("platform_tag", "A")
    label_b = machine_b.get("platform_tag", "B")

    print("=== 环境对比 ===")
    print(f"  文件 A: {path_a}")
    print(f"    机器: {machine_a.get('hostname', '?')}  "
          f"({machine_a.get('platform_tag', '?')})")
    print(f"    快照时间: {snap_a.get('generated_at', '?')}")
    print(f"  文件 B: {path_b}")
    print(f"    机器: {machine_b.get('hostname', '?')}  "
          f"({machine_b.get('platform_tag', '?')})")
    print(f"    快照时间: {snap_b.get('generated_at', '?')}")

    envs_a = snap_a.get("environments", {})
    envs_b = snap_b.get("environments", {})

    total_diffs = 0

    # 对两个环境分别对比
    for env_name in [PY_ENV, R_ENV]:
        total_diffs += _compare_one_env(
            env_name,
            envs_a.get(env_name),
            envs_b.get(env_name),
            label_a,
            label_b,
        )

    # ---- 末尾：如何对齐提示 ----
    print(f"\n{'=' * 70}")
    print("  如何对齐")
    print(f"{'=' * 70}")
    print(
        "  差异需人工决定对齐方向。通常做法：\n"
        "\n"
        "  1. 选择一台机器作为版本基准（通常是已验证可跑的 Mac）。\n"
        "  2. 在落后那台上，对版本不一致的包执行：\n"
        "       conda install -n <env> <pkg>=<基准版本>\n"
        "  3. 对缺失的包同样 `conda install`（或编辑 environment.yml 后重建环境）。\n"
        "  4. 对齐后重新跑 `python scripts/env_parity.py snapshot` 更新快照，\n"
        "     提交到 git 供下次对比。\n"
        "\n"
        "  本脚本不会自动修改环境——对齐决策归人。"
    )

    if total_diffs == 0:
        print("\n✓ 两个环境完全一致，无需对齐。")
        return 0
    else:
        print(f"\n总计 {total_diffs} 项差异，请人工决定对齐方案。")
        return 1


# ============================================================================
# CLI 入口
# ============================================================================


def main() -> int:
    parser = argparse.ArgumentParser(
        description="跨机器环境对齐诊断脚本 —— 导出快照 / 对比差异",
    )
    sub = parser.add_subparsers(dest="command", help="子命令")

    # snapshot
    sub.add_parser("snapshot", help="导出当前机器 conda 环境包清单为 JSON 快照")

    # compare
    cmp_parser = sub.add_parser("compare", help="对比两份快照，列出差异")
    cmp_parser.add_argument(
        "--a",
        default=None,
        help=f"快照文件 A 的路径（默认 {DEFAULT_A}）",
    )
    cmp_parser.add_argument(
        "--b",
        default=None,
        help=f"快照文件 B 的路径（默认 {DEFAULT_B}）",
    )

    args = parser.parse_args()

    if args.command == "snapshot":
        return cmd_snapshot()
    elif args.command == "compare":
        return cmd_compare(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
