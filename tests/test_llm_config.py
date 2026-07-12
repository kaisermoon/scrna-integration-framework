"""Tests for llm_config.py — LLM group config reader.

Uses temporary .env files created via tmp_path fixtures.
No real .env or network calls involved.
"""

from scrna_integration.llm_config import (
    _parse_dotenv,
    _find_vault_root,
    _resolve_endpoint,
    load_llm_group_config,
    get_active_groups,
    call_llm_for_annotation,
    extract_json_from_llm_response,
    build_mllmcelltype_config,
)

import importlib.util

import pytest
from unittest.mock import patch, MagicMock

# call_llm_for_annotation 内部按 provider 惰性 import requests / openai。
# CI 的测试环境未必装这两个库（它们属运行期依赖，非测试依赖），
# 缺库时对应用例无法 patch，跳过而非报错。
_HAS_REQUESTS = importlib.util.find_spec("requests") is not None
_HAS_OPENAI = importlib.util.find_spec("openai") is not None


# ---------------------------------------------------------------------------
# _parse_dotenv
# ---------------------------------------------------------------------------


class TestParseDotenv:
    def test_basic_kv(self, tmp_path):
        p = tmp_path / ".env"
        p.write_text("KEY1=value1\nKEY2=value2\n")
        result = _parse_dotenv(str(p))
        assert result == {"KEY1": "value1", "KEY2": "value2"}

    def test_quoted_value_strips_quotes(self, tmp_path):
        p = tmp_path / ".env"
        p.write_text('KEY="value"\n')
        result = _parse_dotenv(str(p))
        assert result == {"KEY": "value"}

    def test_empty_value_quotes(self, tmp_path):
        p = tmp_path / ".env"
        p.write_text('KEY=""\n')
        result = _parse_dotenv(str(p))
        assert result == {"KEY": ""}

    def test_skips_comments_and_blanks(self, tmp_path):
        p = tmp_path / ".env"
        p.write_text("# comment\n\nKEY=value\n  \n# another\n")
        result = _parse_dotenv(str(p))
        assert result == {"KEY": "value"}

    def test_trims_whitespace(self, tmp_path):
        p = tmp_path / ".env"
        p.write_text("  KEY  =  value  \n")
        result = _parse_dotenv(str(p))
        assert result == {"KEY": "value"}

    def test_multiline_values_not_supported(self, tmp_path):
        """Multi-line values are not supported — only first KEY=VALUE per line."""
        p = tmp_path / ".env"
        p.write_text("KEY=line1\n  line2\n")
        result = _parse_dotenv(str(p))
        # "  line2" has no "=" so it's silently skipped
        assert result == {"KEY": "line1"}

    def test_single_quoted_value_strips_quotes(self, tmp_path):
        """Single-quoted values like KEY='value' should be stripped."""
        p = tmp_path / ".env"
        p.write_text("KEY='value'\n")
        result = _parse_dotenv(str(p))
        assert result == {"KEY": "value"}

    def test_inline_comment_ignored(self, tmp_path):
        """Inline # after space is treated as comment; # inside quotes is preserved."""
        p = tmp_path / ".env"
        p.write_text('KEY1=value # comment\nKEY2="a#b"\n')
        result = _parse_dotenv(str(p))
        assert result == {"KEY1": "value", "KEY2": "a#b"}


# ---------------------------------------------------------------------------
# _find_vault_root
# ---------------------------------------------------------------------------


class TestFindVaultRoot:
    def test_finds_env_in_same_dir(self, tmp_path):
        env = tmp_path / ".env"
        env.write_text("LLM_DEFAULT_GROUP=1\n")
        assert _find_vault_root(str(tmp_path)) == str(tmp_path)

    def test_finds_env_in_parent(self, tmp_path):
        env = tmp_path / ".env"
        env.write_text("LLM_DEFAULT_GROUP=1\n")
        child = tmp_path / "项目" / "my-project"
        child.mkdir(parents=True)
        assert _find_vault_root(str(child)) == str(tmp_path)

    def test_no_env_returns_none(self, tmp_path):
        assert _find_vault_root(str(tmp_path)) is None

    def test_env_without_llm_default_group(self, tmp_path):
        env = tmp_path / ".env"
        env.write_text("OTHER_KEY=value\n")
        assert _find_vault_root(str(tmp_path)) is None


# ---------------------------------------------------------------------------
# load_llm_group_config
# ---------------------------------------------------------------------------


_MINIMAL_ENV = """LLM_DEFAULT_GROUP=1
LLM_GROUP1_PROVIDER=anthropic
LLM_GROUP1_BASE_URL=http://127.0.0.1:8082
LLM_GROUP1_API_KEY=""
LLM_GROUP1_MODEL_HAIKU=claude-haiku-4-5-20251001
LLM_GROUP1_MODEL_SONNET=claude-sonnet-4-6
LLM_GROUP1_MODEL_OPUS=claude-opus-4-8
"""


class TestLoadGroupConfig:
    @staticmethod
    def _make_env(tmp_path, content=_MINIMAL_ENV):
        env = tmp_path / ".env"
        env.write_text(content)
        return str(tmp_path)

    def test_default_group(self, tmp_path):
        root = self._make_env(tmp_path)
        cfg = load_llm_group_config(project_root=root)
        assert cfg is not None
        assert cfg["group"] == 1
        assert cfg["provider"] == "anthropic"
        assert cfg["base_url"] == "http://127.0.0.1:8082"
        assert cfg["api_key"] == ""
        assert cfg["models"]["haiku"] == "claude-haiku-4-5-20251001"
        assert cfg["models"]["sonnet"] == "claude-sonnet-4-6"
        assert cfg["models"]["opus"] == "claude-opus-4-8"
        assert cfg["is_configured"] is True
        assert cfg["has_api_key"] is False  # empty key in _MINIMAL_ENV

    def test_explicit_group_number(self, tmp_path):
        root = self._make_env(tmp_path)
        cfg = load_llm_group_config(group=1, project_root=root)
        assert cfg is not None
        assert cfg["group"] == 1

    def test_unconfigured_group_returns_none(self, tmp_path):
        root = self._make_env(tmp_path)
        cfg = load_llm_group_config(group=2, project_root=root)
        assert cfg is None  # Group2 has empty provider/base_url

    def test_default_group_fallback(self, tmp_path):
        env = tmp_path / ".env"
        env.write_text("LLM_DEFAULT_GROUP=2\n"
                       "LLM_GROUP2_PROVIDER=openai\n"
                       "LLM_GROUP2_BASE_URL=https://api.openai.com\n")
        cfg = load_llm_group_config(project_root=str(tmp_path))
        assert cfg is not None
        assert cfg["group"] == 2
        assert cfg["provider"] == "openai"

    def test_project_root_none_raises(self, tmp_path):
        """Without a .env in the tree, should raise FileNotFoundError."""
        env = tmp_path / ".env"
        env.write_text("NO_LLM_KEY=1\n")  # missing LLM_DEFAULT_GROUP
        with __import__("pytest").raises(FileNotFoundError):
            load_llm_group_config(project_root=str(tmp_path))

    def test_empty_api_key_preserved(self, tmp_path):
        root = self._make_env(tmp_path)
        cfg = load_llm_group_config(project_root=root)
        assert cfg["api_key"] == ""  # empty key for local gateway
        assert cfg["has_api_key"] is False  # empty key = not usable for requests

    def test_api_key_with_value(self, tmp_path):
        env = tmp_path / ".env"
        env.write_text(
            "LLM_DEFAULT_GROUP=1\n"
            "LLM_GROUP1_PROVIDER=openai\n"
            "LLM_GROUP1_BASE_URL=https://api.openai.com\n"
            "LLM_GROUP1_API_KEY=sk-test123\n"
        )
        cfg = load_llm_group_config(project_root=str(tmp_path))
        assert cfg["api_key"] == "sk-test123"
        assert cfg["has_api_key"] is True

    def test_empty_key_detectable_while_group_still_active(self, tmp_path):
        """Bugfix regression: api_key 为空时 group 仍出现在 get_active_groups，
        但 has_api_key=False 让调用方能区分"配置框架在但 key 没填"。
        """
        env = tmp_path / ".env"
        env.write_text(
            "LLM_DEFAULT_GROUP=1\n"
            "LLM_GROUP1_PROVIDER=anthropic\n"
            "LLM_GROUP1_BASE_URL=http://127.0.0.1:8082\n"
            "LLM_GROUP1_API_KEY=\n"  # 空 key —— 缺陷 1 的触发条件
            "LLM_GROUP1_MODEL_HAIKU=claude-haiku-4-5\n"
        )
        # group 仍在 active list 中（provider+base_url 已填）
        active = get_active_groups(project_root=str(tmp_path))
        assert active == [1]

        # 但 has_api_key=False 明确标记 key 缺失
        cfg = load_llm_group_config(group=1, project_root=str(tmp_path))
        assert cfg is not None
        assert cfg["is_configured"] is True
        assert cfg["has_api_key"] is False
        assert cfg["api_key"] == ""

    def test_non_numeric_default_group_fallback(self, tmp_path):
        """Non-numeric LLM_DEFAULT_GROUP falls back to group 1 gracefully."""
        env = tmp_path / ".env"
        env.write_text(
            "LLM_DEFAULT_GROUP=abc\n"
            "LLM_GROUP1_PROVIDER=anthropic\n"
            "LLM_GROUP1_BASE_URL=http://127.0.0.1:8082\n"
        )
        cfg = load_llm_group_config(project_root=str(tmp_path))
        assert cfg is not None
        assert cfg["group"] == 1  # fallback to group 1
        assert cfg["provider"] == "anthropic"

    def test_project_root_none_uses_cwd(self, monkeypatch, tmp_path):
        """project_root=None 时从 os.getcwd() 向上查找 vault 根。

        tmp_path 下无 .env → 应抛出 FileNotFoundError（证明走了 cwd 路径）。
        """
        monkeypatch.setattr("os.getcwd", lambda: str(tmp_path))
        with pytest.raises(FileNotFoundError):
            load_llm_group_config(project_root=None)


class TestGetActiveGroups:
    def test_single_active_group(self, tmp_path):
        env = tmp_path / ".env"
        env.write_text(_MINIMAL_ENV)
        active = get_active_groups(project_root=str(tmp_path))
        assert active == [1]

    def test_multiple_active_groups(self, tmp_path):
        env = tmp_path / ".env"
        env.write_text(
            "LLM_DEFAULT_GROUP=1\n"
            "LLM_GROUP1_PROVIDER=anthropic\n"
            "LLM_GROUP1_BASE_URL=http://127.0.0.1:8082\n"
            "LLM_GROUP2_PROVIDER=openai\n"
            "LLM_GROUP2_BASE_URL=https://api.openai.com\n"
        )
        active = get_active_groups(project_root=str(tmp_path))
        assert active == [1, 2]

    def test_no_active_groups(self, tmp_path):
        env = tmp_path / ".env"
        env.write_text("LLM_DEFAULT_GROUP=1\n"
                       "LLM_GROUP1_PROVIDER=\n"
                       "LLM_GROUP1_BASE_URL=\n")
        active = get_active_groups(project_root=str(tmp_path))
        assert active == []

    def test_project_root_none_uses_cwd(self, monkeypatch, tmp_path):
        """project_root=None 时从 os.getcwd() 向上查找 → 找不到 .env 则抛错。"""
        monkeypatch.setattr("os.getcwd", lambda: str(tmp_path))
        with pytest.raises(FileNotFoundError):
            get_active_groups(project_root=None)


# ---------------------------------------------------------------------------
# extract_json_from_llm_response
# ---------------------------------------------------------------------------


class TestExtractJsonFromLLMResponse:
    def test_pure_json_dict(self):
        result = extract_json_from_llm_response('{"cell_type": "B cell", "score": "0.95"}')
        assert result == {"cell_type": "B cell", "score": "0.95"}

    def test_json_code_fence(self):
        raw = '```json\n{"cell_type": "T cell"}\n```'
        result = extract_json_from_llm_response(raw)
        assert result == {"cell_type": "T cell"}

    def test_code_fence_without_lang(self):
        raw = '```\n{"cell_type": "NK cell"}\n```'
        result = extract_json_from_llm_response(raw)
        assert result == {"cell_type": "NK cell"}

    def test_code_fence_with_prefix_text(self):
        raw = 'Here is the annotation:\n```json\n{"cell_type": "B cell"}\n```\nHope this helps.'
        result = extract_json_from_llm_response(raw)
        assert result == {"cell_type": "B cell"}

    def test_embedded_dict_no_fence(self):
        raw = 'The result is {"cell_type": "Macrophage"} for this cluster.'
        result = extract_json_from_llm_response(raw)
        assert result == {"cell_type": "Macrophage"}

    def test_empty_input_returns_none(self):
        assert extract_json_from_llm_response("") is None

    def test_non_json_text_returns_none(self):
        assert extract_json_from_llm_response("This is not JSON at all.") is None

    def test_array_not_dict_returns_none(self):
        raw = '```json\n["a", "b", "c"]\n```'
        result = extract_json_from_llm_response(raw)
        assert result is None

    def test_filters_empty_and_null_values(self):
        raw = '{"cell_type": "B cell", "empty": "", "null_val": null, "score": "0.9"}'
        result = extract_json_from_llm_response(raw)
        assert result == {"cell_type": "B cell", "score": "0.9"}

    def test_nested_json_is_valid_dict(self):
        """嵌套 dict 是合法 JSON，json.loads 正常解析后应成功返回 dict。"""
        raw = '{"outer": {"inner": "value"}}'
        result = extract_json_from_llm_response(raw)
        # 嵌套 dict 是合法 dict，应成功返回
        assert isinstance(result, dict)
        assert "outer" in result


# ---------------------------------------------------------------------------
# _resolve_endpoint
# ---------------------------------------------------------------------------


class TestResolveEndpoint:
    def test_v1_slash_already_present(self):
        assert _resolve_endpoint("http://127.0.0.1:8082/v1/", "anthropic") == "http://127.0.0.1:8082/v1/"

    def test_empty_url(self):
        assert _resolve_endpoint("", "anthropic") == ""

    def test_unknown_provider_no_suffix(self):
        assert _resolve_endpoint("http://127.0.0.1:8082", "unknown_provider") == "http://127.0.0.1:8082"

    def test_trailing_slash_handled(self):
        assert _resolve_endpoint("http://127.0.0.1:8082/", "anthropic") == "http://127.0.0.1:8082/v1/messages"

    def test_anthropic_endpoint(self):
        assert _resolve_endpoint("http://127.0.0.1:8082", "anthropic") == "http://127.0.0.1:8082/v1/messages"

    def test_openai_endpoint(self):
        assert _resolve_endpoint("https://api.openai.com", "openai") == "https://api.openai.com/v1"

    def test_v1_no_trailing_slash_no_double_write(self):
        """Issue #7 regression: base_url="/v1" (无尾斜杠) 不应双写成 /v1/v1/messages。"""
        assert _resolve_endpoint("http://127.0.0.1:8082/v1", "anthropic") == "http://127.0.0.1:8082/v1"

    def test_v1_messages_already_complete(self):
        assert _resolve_endpoint("http://127.0.0.1:8082/v1/messages", "anthropic") == "http://127.0.0.1:8082/v1/messages"


# ---------------------------------------------------------------------------
# call_llm_for_annotation
# ---------------------------------------------------------------------------


_ENV_ANTHROPIC = """LLM_DEFAULT_GROUP=1
LLM_GROUP1_PROVIDER=anthropic
LLM_GROUP1_BASE_URL=http://127.0.0.1:8082
LLM_GROUP1_API_KEY=test-key
LLM_GROUP1_MODEL_HAIKU=claude-haiku-4-5
"""

_ENV_OPENAI = """LLM_DEFAULT_GROUP=1
LLM_GROUP1_PROVIDER=openai
LLM_GROUP1_BASE_URL=https://api.openai.com
LLM_GROUP1_API_KEY=sk-test
LLM_GROUP1_MODEL_HAIKU=gpt-4o-mini
"""


@pytest.mark.skipif(
    not (_HAS_REQUESTS and _HAS_OPENAI),
    reason="call_llm_for_annotation 测试需 requests 和 openai（运行期依赖，CI 环境可选）",
)
class TestCallLLMForAnnotation:
    @staticmethod
    def _make_env(tmp_path, content):
        env = tmp_path / ".env"
        env.write_text(content)
        return str(tmp_path)

    def test_anthropic_happy_path(self, tmp_path):
        root = self._make_env(tmp_path, _ENV_ANTHROPIC)
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "content": [{"type": "text", "text": "B cell"}, {"type": "text", "text": "T cell"}]
        }
        with patch("requests.post", return_value=mock_resp) as mock_post:
            result = call_llm_for_annotation(
                [{"role": "user", "content": "Annotate"}],
                project_root=root,
            )
        assert result == "B cell\nT cell"
        mock_post.assert_called_once()
        # 验证 raise_for_status 被调用（issue #5 修复）
        mock_resp.raise_for_status.assert_called_once()

    def test_anthropic_error_raises(self, tmp_path):
        root = self._make_env(tmp_path, _ENV_ANTHROPIC)
        mock_resp = MagicMock()
        mock_resp.raise_for_status.side_effect = __import__("requests").exceptions.HTTPError("401 Unauthorized")
        with patch("requests.post", return_value=mock_resp):
            with pytest.raises(__import__("requests").exceptions.HTTPError, match="401"):
                call_llm_for_annotation(
                    [{"role": "user", "content": "Annotate"}],
                    project_root=root,
                )

    def test_openai_happy_path(self, tmp_path):
        root = self._make_env(tmp_path, _ENV_OPENAI)
        mock_client = MagicMock()
        mock_choice = MagicMock()
        mock_choice.message.content = "B cell"
        mock_client.chat.completions.create.return_value.choices = [mock_choice]
        with patch("openai.OpenAI", return_value=mock_client):
            result = call_llm_for_annotation(
                [{"role": "user", "content": "Annotate"}],
                project_root=root,
            )
        assert result == "B cell"

    def test_openai_error_propagates(self, tmp_path):
        root = self._make_env(tmp_path, _ENV_OPENAI)
        mock_client = MagicMock()
        mock_client.chat.completions.create.side_effect = __import__("openai").APIError(
            "Rate limit", request=None, body=None
        )
        with patch("openai.OpenAI", return_value=mock_client):
            with pytest.raises(__import__("openai").APIError):
                call_llm_for_annotation(
                    [{"role": "user", "content": "Annotate"}],
                    project_root=root,
                )

    def test_config_provider_overrides_model_provider_param(self, tmp_path):
        """当 group 配置已指定 provider（anthropic）时，忽略 model_provider 参数（openai）。"""
        root = self._make_env(tmp_path, _ENV_ANTHROPIC)
        mock_resp = MagicMock()
        mock_resp.json.return_value = {"content": [{"type": "text", "text": "OK"}]}
        with patch("requests.post", return_value=mock_resp):
            result = call_llm_for_annotation(
                [{"role": "user", "content": "Annotate"}],
                model_provider="openai",  # fallback — config has "anthropic" so this is ignored
                project_root=root,
            )
        assert result == "OK"


# ---------------------------------------------------------------------------
# build_mllmcelltype_config
# ---------------------------------------------------------------------------


class TestBuildMLLMCelltypeConfig:
    @staticmethod
    def _make_env(tmp_path, content):
        env = tmp_path / ".env"
        env.write_text(content)
        return str(tmp_path)

    def test_single_group_triple(self, tmp_path):
        root = self._make_env(tmp_path, _ENV_ANTHROPIC)
        api_keys, base_urls, model_list = build_mllmcelltype_config(project_root=root)
        assert api_keys == {"anthropic": "test-key"}
        assert base_urls == {"anthropic": "http://127.0.0.1:8082/v1/messages"}
        assert model_list == ["claude-haiku-4-5"]

    def test_multiple_groups(self, tmp_path):
        env = tmp_path / ".env"
        env.write_text(
            "LLM_DEFAULT_GROUP=1\n"
            "LLM_GROUP1_PROVIDER=anthropic\n"
            "LLM_GROUP1_BASE_URL=http://127.0.0.1:8082\n"
            "LLM_GROUP1_API_KEY=key1\n"
            "LLM_GROUP1_MODEL_HAIKU=claude-haiku-4-5\n"
            "LLM_GROUP2_PROVIDER=openai\n"
            "LLM_GROUP2_BASE_URL=https://api.openai.com\n"
            "LLM_GROUP2_API_KEY=key2\n"
            "LLM_GROUP2_MODEL_HAIKU=gpt-4o-mini\n"
        )
        api_keys, base_urls, model_list = build_mllmcelltype_config(project_root=str(tmp_path))
        assert api_keys == {"anthropic": "key1", "openai": "key2"}
        assert len(base_urls) == 2
        assert set(model_list) == {"claude-haiku-4-5", "gpt-4o-mini"}

    def test_model_list_override_with_keys(self, tmp_path):
        """Issue #6 修复：有 api_key 时 override 应被接受，不应误报 missing。"""
        root = self._make_env(tmp_path, _ENV_ANTHROPIC)
        api_keys, base_urls, model_list = build_mllmcelltype_config(
            project_root=root,
            model_list_override=["claude-opus-4-8", "gpt-4o"],
        )
        # api_keys 非空，override 应被完整接受
        assert model_list == ["claude-opus-4-8", "gpt-4o"]
        assert api_keys == {"anthropic": "test-key"}

    def test_model_list_override_without_keys(self, tmp_path):
        """Issue #6 修复：无 api_key 时 override 应全部报告 missing 并返回空列表。"""
        env = tmp_path / ".env"
        env.write_text(
            "LLM_DEFAULT_GROUP=1\n"
            "LLM_GROUP1_PROVIDER=anthropic\n"
            "LLM_GROUP1_BASE_URL=http://127.0.0.1:8082\n"
            "LLM_GROUP1_API_KEY=\n"  # 空 key
            "LLM_GROUP1_MODEL_HAIKU=claude-haiku-4-5\n"
        )
        api_keys, base_urls, model_list = build_mllmcelltype_config(
            project_root=str(tmp_path),
            model_list_override=["claude-opus-4-8"],
        )
        # 空 api_key 不会被加入 api_keys → api_keys 为空 → override 全部 missing
        assert model_list == []
        # api_keys 应为空（key 为空字符串时不会加入）
        assert api_keys == {}

    def test_no_active_groups(self, tmp_path):
        env = tmp_path / ".env"
        env.write_text("LLM_DEFAULT_GROUP=1\n"
                       "LLM_GROUP1_PROVIDER=\n"
                       "LLM_GROUP1_BASE_URL=\n")
        api_keys, base_urls, model_list = build_mllmcelltype_config(project_root=str(tmp_path))
        assert api_keys == {}
        assert base_urls == {}
        assert model_list == []

    def test_model_list_override_with_no_active_groups(self, tmp_path):
        """Issue #6 修复：无活跃 group 时 override 应全部 missing。"""
        env = tmp_path / ".env"
        env.write_text("LLM_DEFAULT_GROUP=1\n"
                       "LLM_GROUP1_PROVIDER=\n"
                       "LLM_GROUP1_BASE_URL=\n")
        api_keys, base_urls, model_list = build_mllmcelltype_config(
            project_root=str(tmp_path),
            model_list_override=["claude-opus-4-8"],
        )
        assert model_list == []

    def test_b7_multi_tier_single_group_all_three(self, tmp_path):
        """B7 行为变更：单个 group 配满 haiku + sonnet + opus → model_list 含 3 个模型。

        这验证 build_mllmcelltype_config 从"每 group 取 haiku 单模型"
        改为"取每个 group 所有已配置档位"的新行为。
        """
        env = tmp_path / ".env"
        env.write_text(
            "LLM_DEFAULT_GROUP=1\n"
            "LLM_GROUP1_PROVIDER=anthropic\n"
            "LLM_GROUP1_BASE_URL=http://127.0.0.1:8082\n"
            "LLM_GROUP1_API_KEY=test-key\n"
            "LLM_GROUP1_MODEL_HAIKU=claude-haiku-4-5\n"
            "LLM_GROUP1_MODEL_SONNET=claude-sonnet-4-6\n"
            "LLM_GROUP1_MODEL_OPUS=claude-opus-4-8\n"
        )
        api_keys, base_urls, model_list = build_mllmcelltype_config(
            project_root=str(tmp_path),
        )
        assert len(model_list) == 3
        assert model_list == [
            "claude-haiku-4-5",
            "claude-sonnet-4-6",
            "claude-opus-4-8",
        ]
        assert api_keys == {"anthropic": "test-key"}
        assert "anthropic" in base_urls

    def test_b7_multi_tier_single_group_haiku_only(self, tmp_path):
        """B7 行为变更：单 group 只配 haiku → model_list 仅含 1 个模型（兼容旧行为）。

        确保没有 sonnet/opus 时不会引入空字符串或假模型名。
        """
        root = self._make_env(tmp_path, _ENV_ANTHROPIC)
        api_keys, base_urls, model_list = build_mllmcelltype_config(
            project_root=root,
        )
        assert len(model_list) == 1
        assert model_list == ["claude-haiku-4-5"]
        assert api_keys == {"anthropic": "test-key"}
        assert "anthropic" in base_urls
