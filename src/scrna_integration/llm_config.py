"""LLM group 配置读取。

从 vault 根 .env 的 ``LLM_GROUP{N}_*`` schema 读取 LLM 配置，
替代旧的 ``{PROVIDER}_API_KEY`` 散落写法。
一个 group = 一个 LLM 端点 + provider + 模型系列。

面向非计算机专业 PI/学生：打开文件从上到下可读。

公开函数：
- ``load_llm_group_config()`` — 读取单个 group 的完整配置
- ``get_active_groups()`` — 列出所有已配置的 group 编号
"""

from __future__ import annotations

import os
from pathlib import Path


# ---------------------------------------------------------------------------
# 内部 helper：vault 根发现 + 手动 .env 解析（零额外依赖）
# ---------------------------------------------------------------------------


def _find_vault_root(project_root: str | Path) -> str | None:
    """从项目根向上查找 AI-OS vault 根。

    vault 根的特征：包含 .env 文件且其中有 ``LLM_DEFAULT_GROUP=`` 行。
    最多向上遍历到文件系统根。
    """
    path = Path(project_root).resolve()
    for parent in [path] + list(path.parents):
        env_file = parent / ".env"
        if not env_file.is_file():
            continue
        try:
            content = env_file.read_text()
            if "LLM_DEFAULT_GROUP=" in content:
                return str(parent)
        except Exception:
            continue
    return None


def _parse_dotenv(env_path: str | Path) -> dict[str, str]:
    """手动解析 .env 为 dict——零依赖，避免引入 python-dotenv。

    规则：
    - 跳过空行和 ``#`` 开头的注释行
    - 解析 ``KEY=VALUE``，去除两端空格
    - VALUE 两端引号（单引号或双引号）自动剥除
    - 行内 `` #`` 注释剥离（引号内的 ``#`` 受保护）

    处理顺序：先判断值是否被引号整体包裹——若被包裹则剥引号后
    ``#`` 属于值的一部分，不剥离；若未被包裹则先剥行内注释，再检查
    剩余部分是否被引号包裹。这样 ``KEY="a#b"`` → ``a#b`` 不会误伤。
    """
    result: dict[str, str] = {}
    with open(env_path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip()
            # 若值整体被引号包裹（双引号或单引号），剥引号后 # 是值的一部分
            if len(value) >= 2 and value[0] == value[-1] and value[0] in ('"', "'"):
                value = value[1:-1]
            else:
                # 未被引号整体包裹 → 行内 # 视为注释起始
                # 用 " #"（空格+#）分隔，避免误伤值内的 # 字符
                value = value.split(" #")[0].rstrip()
                # 剥注释后残留空格可能使值重新被引号包裹（如 KEY='val' # cmt）
                if len(value) >= 2 and value[0] == value[-1] and value[0] in ('"', "'"):
                    value = value[1:-1]
            result[key] = value
    return result


# ---------------------------------------------------------------------------
# 公开函数
# ---------------------------------------------------------------------------


def load_llm_group_config(
    group: int | None = None,
    project_root: str | None = None,
) -> dict | None:
    """读取指定 LLM group 的完整配置。

    Parameters
    ----------
    group : int or None
        要读取的 group 编号（1 起）。None 时自动取 ``LLM_DEFAULT_GROUP``。
    project_root : str or None
        项目根目录。None 时从 ``os.getcwd()`` 向上推断。

    Returns
    -------
    dict or None
        若 group 已配置（provider + base_url 均非空）则返回 dict::

            {
                "group": 1,
                "name": "anthropic",
                "provider": "anthropic",
                "base_url": "http://127.0.0.1:8082",
                "api_key": "",
                "models": {"haiku": "claude-haiku-4-5-20251001",
                           "sonnet": "claude-sonnet-4-6",
                           "opus": "claude-opus-4-8"},
                "is_configured": True,
                "has_api_key": False,
            }

        has_api_key 仅表示配置中有非空 key 字符串，不保证 key 能通过远端认证。

        若 group 未配置（provider/base_url 为空）则返回 None。
    """
    if project_root is None:
        project_root = os.getcwd()

    vault_root = _find_vault_root(project_root)
    if vault_root is None:
        raise FileNotFoundError(
            "找不到 AI-OS vault 根 .env（需包含 LLM_DEFAULT_GROUP= 行）。"
            "请确认项目位于 vault 内正确的位置。"
        )

    env = _parse_dotenv(Path(vault_root) / ".env")

    # 确定 group 编号
    if group is None:
        default_str = env.get("LLM_DEFAULT_GROUP", "1").strip()
        try:
            group = int(default_str)
        except ValueError:
            group = 1

    prefix = f"LLM_GROUP{group}_"

    provider = env.get(f"{prefix}PROVIDER", "").strip()
    base_url = env.get(f"{prefix}BASE_URL", "").strip()
    api_key = env.get(f"{prefix}API_KEY", "").strip()
    name = env.get(f"{prefix}NAME", "").strip()

    models = {
        "haiku": env.get(f"{prefix}MODEL_HAIKU", "").strip(),
        "sonnet": env.get(f"{prefix}MODEL_SONNET", "").strip(),
        "opus": env.get(f"{prefix}MODEL_OPUS", "").strip(),
    }

    is_configured = bool(provider and base_url)
    has_api_key = bool(api_key)

    if not is_configured:
        return None

    return {
        "group": group,
        "name": name or provider,
        "provider": provider,
        "base_url": base_url,
        "api_key": api_key,
        "models": models,
        "is_configured": True,
        "has_api_key": has_api_key,
    }


def get_active_groups(project_root: str | None = None) -> list[int]:
    """返回所有已配置的 group 编号（provider + base_url 均非空）。

    按编号升序。最多检查 1-9 号 group。
    用于 notebook 判断当前有几个 LLM 端点可用。
    """
    active: list[int] = []
    for g in range(1, 10):
        cfg = load_llm_group_config(group=g, project_root=project_root)
        if cfg is not None:
            active.append(g)
    return active
