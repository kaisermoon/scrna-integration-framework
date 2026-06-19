"""Tests for llm_config.py — LLM group config reader.

Uses temporary .env files created via tmp_path fixtures.
No real .env or network calls involved.
"""

from scrna_integration.llm_config import (
    _parse_dotenv,
    _find_vault_root,
    load_llm_group_config,
    get_active_groups,
)

import pytest


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
