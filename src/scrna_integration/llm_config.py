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


# ---------------------------------------------------------------------------
# 公开函数：LLM 统一调用（批 5：收编 06/06b/06c 三处漂移的 LLM 调用逻辑）
# ---------------------------------------------------------------------------

# 按 provider 补端点的映射表
_PROVIDER_ENDPOINT: dict[str, str] = {
    "anthropic": "/v1/messages",
    "openai": "/v1",
    "deepseek": "/v1",
    "groq": "/v1",
    "together": "/v1",
    "fireworks": "/v1",
    "siliconflow": "/v1",
}


def _resolve_endpoint(base_url: str, provider: str) -> str:
    """补全 base_url 的 API 端点后缀。

    若 base_url 已含 ``/v1/`` 路径或以 ``/v1`` 结尾则不补，
    避免双写（如 ``/v1/v1/messages``）。
    否则按 provider 查 ``_PROVIDER_ENDPOINT`` 追加对应后缀。
    """
    _normalized = base_url.rstrip("/")
    if not base_url or "/v1/" in base_url or _normalized.endswith("/v1"):
        return base_url
    suffix = _PROVIDER_ENDPOINT.get(provider, "")
    if not suffix:
        return base_url
    base = base_url.rstrip("/")
    return f"{base}{suffix}"


def call_llm_for_annotation(
    messages: list[dict],
    *,
    model_provider: str = "anthropic",
    model_name: str = "",
    max_tokens: int = 1024,
    temperature: float = 0.2,
    system: str = "",
    project_root: str | None = None,
    group: int | None = None,
    timeout: int = 120,
) -> str:
    """统一的 LLM 调用函数。

    从 vault 根 .env 读取 LLM 配置，按 provider 分支调用对应 SDK，
    返回响应中的文本内容。

    Parameters
    ----------
    messages : list[dict]
        消息列表，格式为 ``[{"role": "user", "content": "..."}]``。
    model_provider : str
        ``"anthropic"`` 或 ``"openai"``。
        仅当 group 配置中未指定 provider 时作为 fallback 生效；
        若配置中已设 provider，以配置为准，此参数被忽略。
    model_name : str
        具体模型名（如 ``"claude-haiku-4-5-20251001"``）。
        为空时自动取 group 配置中 haiku 档模型。
    max_tokens : int
        最大输出 token 数。
    temperature : float
        采样温度（0-1）。
    system : str
        系统提示词。anthropic 作为顶层 system 字段发送；
        openai-compatible 作为 messages 中的 system role 消息。
    project_root : str or None
        项目根目录。None 时自动推断。
    group : int or None
        LLM group 编号。None 时取 LLM_DEFAULT_GROUP。
    timeout : int
        HTTP 请求超时秒数。

    Returns
    -------
    str
        LLM 返回的文本内容。
    """
    cfg = load_llm_group_config(group=group, project_root=project_root)
    if cfg is None:
        raise RuntimeError(
            "LLM group 未配置。请在 vault 根 .env 中设置 "
            "LLM_GROUP{N}_PROVIDER / LLM_GROUP{N}_BASE_URL。"
        )

    api_key: str = cfg.get("api_key", "")
    base_url: str = cfg.get("base_url", "")
    provider: str = cfg.get("provider", model_provider)
    models: dict = cfg.get("models", {})

    if not model_name:
        model_name = models.get("haiku", "")
    if not model_name:
        raise ValueError(
            "model_name 未指定，且 group 配置中无 haiku 档模型。"
        )

    # 端点补全
    resolved_url = _resolve_endpoint(base_url, provider)

    # 按 provider 分支调用
    if provider == "anthropic":
        import requests as _req

        _headers = {
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01",
            "content-type": "application/json",
        }
        _body: dict = {
            "model": model_name,
            "max_tokens": max_tokens,
            "messages": messages,
        }
        if system:
            _body["system"] = system

        _resp = _req.post(
            resolved_url,
            headers=_headers,
            json=_body,
            timeout=timeout,
        )
        _resp.raise_for_status()
        _data = _resp.json()
        _text_blocks = [
            b.get("text", "")
            for b in _data.get("content", [])
            if b.get("type") == "text"
        ]
        return "\n".join(_text_blocks)

    else:
        # OpenAI-compatible providers（openai / deepseek / groq / ...）
        from openai import OpenAI as _OpenAI

        _client = _OpenAI(
            api_key=api_key,
            base_url=resolved_url,
            timeout=timeout,
        )
        _api_messages: list[dict] = []
        if system:
            _api_messages.append({"role": "system", "content": system})
        _api_messages.extend(messages)

        _completion = _client.chat.completions.create(
            model=model_name,
            max_tokens=max_tokens,
            temperature=temperature,
            messages=_api_messages,
        )
        return _completion.choices[0].message.content or ""


def extract_json_from_llm_response(raw_text: str) -> dict | None:
    """从 LLM 文本回复中提取 JSON 字典。

    容忍多种格式：纯 JSON、`` ```json ````` 代码围栏、
    嵌入文本中的 JSON 等。

    Parameters
    ----------
    raw_text : str
        LLM 返回的原始文本（可能含前言/后缀）。

    Returns
    -------
    dict or None
        解析成功返回 dict；解析失败返回 None。
    """
    import json as _json
    import re as _re

    if not raw_text:
        return None

    # 策略 1：匹配 ```json ... ``` 代码围栏
    _match = _re.search(
        r'```(?:json)?\s*\n?(.*?)\n?```', raw_text, _re.DOTALL
    )
    if _match:
        _candidate = _match.group(1).strip()
        try:
            _result = _json.loads(_candidate)
            if isinstance(_result, dict):
                return {str(k): str(v) for k, v in _result.items() if v}
        except (_json.JSONDecodeError, ValueError):
            pass

    # 策略 2：找第一个 { 和最后一个 } 之间的内容
    _start = raw_text.find("{")
    _end = raw_text.rfind("}")
    if _start >= 0 and _end > _start:
        _candidate = raw_text[_start:_end + 1]
        try:
            _result = _json.loads(_candidate)
            if isinstance(_result, dict):
                return {str(k): str(v) for k, v in _result.items() if v}
        except (_json.JSONDecodeError, ValueError):
            pass

    return None


def build_mllmcelltype_config(
    project_root: str | None = None,
    model_list_override: list[str] | None = None,
) -> tuple[dict[str, str], dict[str, str], list[str]]:
    """从 .env LLM_GROUP 构建 mLLMCelltype 所需的配置。

    遍历所有活跃 group，取每个 group 的所有已配置模型档位 +
    api_key + base_url（已补端点）。返回的字典按唯一 group 键组织
    （NAME 优先、未设时回落 provider、同 provider 多个 group 冲突时
    追加 ``#编号`` 去重），不再按 provider 归并。

    Parameters
    ----------
    project_root : str or None
        项目根目录。
    model_list_override : list[str] or None
        手动指定的模型列表（覆盖自动构建）。
        若所有 provider 均无 api_key，则打印警告并返回空列表；
        否则以 override 列表为准（下游库负责模型-提供者路由匹配）。

    Returns
    -------
    tuple
        ``(api_keys, base_urls, model_list)`` 三元组。
        - ``api_keys``: ``{group_key: api_key}``
        - ``base_urls``: ``{group_key: resolved_url}``
        - ``model_list``: 模型名列表

        group_key 构成规则：优先取 group 的 NAME（LLM_GROUP{N}_NAME），
        未设时回落为 provider 名；若该键已被前一个贡献了 api_key 的
        group 占用，追加 ``#编号`` 后缀以保证每个 group 拥有独立键。
    """
    active_groups = get_active_groups(project_root=project_root)
    if not active_groups:
        return {}, {}, []

    api_keys: dict[str, str] = {}
    base_urls: dict[str, str] = {}
    model_list: list[str] = []

    # 记录已被前面 group 占用的路由键，用于同 provider 多个 group
    # 冲突时自动追加 gid 编号去重（决策7）。
    _seen_group_keys: set[str] = set()

    for gid in active_groups:
        cfg = load_llm_group_config(group=gid, project_root=project_root)
        if cfg is None:
            continue
        provider = cfg.get("provider", "")
        key = cfg.get("api_key", "")
        url = cfg.get("base_url", "")

        # 决策7：路由键改用唯一 group 标识，避免多个 group 共用同一
        # provider 时后组覆盖前组的 key 与 endpoint。优先用 group
        # 名（LLM_GROUP{N}_NAME，未设时回落 provider），若该名已被
        # 前一个 group 占用，追加 group 编号去重，保证每个 group
        # 拥有独立键。
        group_key: str = cfg.get("name", "") or provider
        if group_key in _seen_group_keys:
            group_key = f"{group_key}#{gid}"

        # B7 默认多模型共识：取每个 group 所有已配置的模型档位
        # （haiku / sonnet / opus），而非仅 haiku 单模型。
        # 下游调用方根据模型数量判断是否启用 discussion 模式。
        _contributed = False
        for _tier in ("haiku", "sonnet", "opus"):
            _model = cfg.get("models", {}).get(_tier, "")
            if _model and key:
                model_list.append(_model)
                api_keys[group_key] = key
                if url:
                    base_urls[group_key] = _resolve_endpoint(url, provider)
                _contributed = True
        if _contributed:
            _seen_group_keys.add(group_key)

    if model_list_override:
        # api_keys 的键是唯一 group 键（NAME 优先 / provider 回落 /
        # 冲突追加 #编号），而 model_list_override 中的是模型名
        # （如 "claude-haiku-4-5"），不能直接用模型名去匹配 group 键。
        # 若没有任何已配置的 api_key，则所有 override 模型都缺少 key；
        # 否则接受 override 列表，下游路由会将模型分派到对应 provider。
        if not api_keys:
            print(
                f"MLLM_MODELS 中以下模型无 API key: {model_list_override}"
            )
            model_list = []
        else:
            model_list = model_list_override

    return api_keys, base_urls, model_list


def apply_mllmcelltype_patches(
    max_retries: int = 3,
    retry_delay: int = 2,
    timeout: int = 120,
    request_json: bool = False,
    max_tokens_override: int = 16384,
) -> bool:
    """对 mLLMCelltype 库内部函数应用 monkey-patch。

    修复库的已知缺陷：超时过小、max_tokens 过低、
    Anthropic thinking 块解析不兼容。

    参数全部显式命名，面向非 CS 学生可读。

    Returns
    -------
    bool
        True 表示 patch 成功；False 表示库版本不兼容，跳过。
    """
    try:
        import mllmcelltype.functions as _mf
        import mllmcelltype.providers.common as _c
        import mllmcelltype.providers.anthropic as _a
        import requests as _r

        # 修复 1：HTTP 重试默认参数（用显式变量名替代位置元组）
        _c.call_http_api_with_retry.__defaults__ = (
            max_retries,    # max_retries
            retry_delay,    # retry_delay (seconds)
            timeout,        # timeout (seconds)
            request_json,   # request_json
            (),             # extra (empty tuple)
        )

        # 修复 2：OpenAI 兼容 API 默认参数
        _c.call_openai_compatible_api.__defaults__ = (
            max_tokens_override,  # max_tokens
            None,                 # top_p
            (),                   # extra
            None,                 # temperature
            None,                 # top_k
            None,                 # repetition_penalty
            False,                # request_json
        )

        # 修复 3：Anthropic thinking 块解析兼容
        _a._parse_anthropic_response = lambda content: [
            line.rstrip(",")
            for block in content.get("content", [])
            if block.get("type") == "text" and "text" in block
            for line in block["text"].strip().split("\n")
        ]

        # 修复 4：注册 patched Anthropic provider 函数
        def _patched_process_anthropic(prompt, model, api_key,
                                        base_url=None):
            return _c.call_http_api_with_retry(
                provider_name="Anthropic",
                url=_a.resolve_endpoint_url(
                    "anthropic", "Anthropic", base_url
                ),
                body={
                    "model": _a._resolve_model_name(model),
                    "messages": [{"role": "user", "content": prompt}],
                    "max_tokens": max_tokens_override,
                    "thinking": {"type": "disabled"},
                },
                headers={
                    "Content-Type": "application/json",
                    "x-api-key": api_key,
                    "anthropic-version": "2023-06-01",
                },
                post_func=_r.post,
                response_parser=_a._parse_anthropic_response,
                max_retries=max_retries,
                retry_delay=retry_delay,
                timeout=timeout,
                request_json=False,
            )

        _mf.PROVIDER_FUNCTIONS["anthropic"] = (
            _patched_process_anthropic
        )
        return True

    except (AttributeError, ImportError, TypeError) as e:
        import warnings
        warnings.warn(
            f"mLLMCelltype monkey-patch 失败（库版本可能不兼容）: {e}\n"
            "  mLLMCelltype 注释将跳过。"
            "  若需恢复，请检查 mllmcelltype 版本"
            "（当前验证版本: 2.1.1）。"
        )
        return False
