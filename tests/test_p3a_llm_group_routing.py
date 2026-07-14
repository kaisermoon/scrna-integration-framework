"""P3-a 叶子测试：build_mllmcelltype_config() 路由键从 provider 改为唯一 group 标识。

修复多 group 同 provider 时后组覆盖前组 api_key/base_url 的 bug（决策7）。
全部用例用 tmp_path 造临时 .env，零网络依赖。
"""

from scrna_integration.llm_config import build_mllmcelltype_config


class TestGroupRoutingKey:
    """验证 build_mllmcelltype_config 返回的 api_keys/base_urls 键语义。"""

    @staticmethod
    def _make_env(tmp_path, content: str) -> str:
        env = tmp_path / ".env"
        env.write_text(content)
        return str(tmp_path)

    # ------------------------------------------------------------------
    # 核心：同 provider 多 group 不再互相覆盖
    # ------------------------------------------------------------------

    def test_same_provider_distinct_names_both_preserved(self, tmp_path):
        """两个 group 同 provider=anthropic，各自设 NAME → 两键均保留，不覆盖。"""
        env_content = """LLM_DEFAULT_GROUP=1
LLM_GROUP1_PROVIDER=anthropic
LLM_GROUP1_BASE_URL=http://127.0.0.1:8082
LLM_GROUP1_API_KEY=key-local
LLM_GROUP1_NAME=local-gateway
LLM_GROUP1_MODEL_HAIKU=claude-haiku-4-5
LLM_GROUP2_PROVIDER=anthropic
LLM_GROUP2_BASE_URL=https://api.anthropic.com
LLM_GROUP2_API_KEY=key-cloud
LLM_GROUP2_NAME=cloud-direct
LLM_GROUP2_MODEL_HAIKU=claude-haiku-4-5
"""
        root = self._make_env(tmp_path, env_content)
        api_keys, base_urls, model_list = build_mllmcelltype_config(project_root=root)

        # 键为各自 NAME，两条都保留
        assert api_keys == {"local-gateway": "key-local", "cloud-direct": "key-cloud"}
        assert len(base_urls) == 2
        assert "local-gateway" in base_urls
        assert "cloud-direct" in base_urls
        # 各自 url 均被保留（互不覆盖）
        assert base_urls["local-gateway"].startswith("http://127.0.0.1:8082")
        assert base_urls["cloud-direct"].startswith("https://api.anthropic.com")
        assert len(model_list) == 2

    def test_same_provider_no_names_disambiguated(self, tmp_path):
        """两个 group 同 anthropic、均不设 NAME → 一个键为 provider 名，另一个追加 #2。"""
        env_content = """LLM_DEFAULT_GROUP=1
LLM_GROUP1_PROVIDER=anthropic
LLM_GROUP1_BASE_URL=http://127.0.0.1:8082
LLM_GROUP1_API_KEY=key1
LLM_GROUP1_MODEL_HAIKU=claude-haiku-4-5
LLM_GROUP2_PROVIDER=anthropic
LLM_GROUP2_BASE_URL=http://127.0.0.1:8083
LLM_GROUP2_API_KEY=key2
LLM_GROUP2_MODEL_HAIKU=claude-haiku-4-5
"""
        root = self._make_env(tmp_path, env_content)
        api_keys, base_urls, model_list = build_mllmcelltype_config(project_root=root)

        # 两条 key 都保留（值不丢）
        assert len(api_keys) == 2
        assert set(api_keys.values()) == {"key1", "key2"}
        # 键：一个为裸 "anthropic"，另一个含 "#" 去重后缀
        keys = list(api_keys.keys())
        assert "anthropic" in keys
        assert any(k.startswith("anthropic#") for k in keys)

        # base_urls 同样两条都在
        assert len(base_urls) == 2
        assert set(base_urls.keys()) == set(api_keys.keys())

    # ------------------------------------------------------------------
    # 回归锁：单 provider / 异 provider 保持与现状逐字节一致
    # ------------------------------------------------------------------

    def test_single_provider_key_unchanged(self, tmp_path):
        """单个 anthropic group → 键仍为 provider 名 'anthropic'。"""
        env_content = """LLM_DEFAULT_GROUP=1
LLM_GROUP1_PROVIDER=anthropic
LLM_GROUP1_BASE_URL=http://127.0.0.1:8082
LLM_GROUP1_API_KEY=test-key
LLM_GROUP1_MODEL_HAIKU=claude-haiku-4-5
"""
        root = self._make_env(tmp_path, env_content)
        api_keys, base_urls, model_list = build_mllmcelltype_config(project_root=root)

        assert api_keys == {"anthropic": "test-key"}
        assert base_urls == {"anthropic": "http://127.0.0.1:8082/v1/messages"}
        assert model_list == ["claude-haiku-4-5"]

    def test_distinct_providers_unchanged(self, tmp_path):
        """anthropic + openai 两 group → 键仍为 provider 名，无 '#' 后缀。"""
        env_content = """LLM_DEFAULT_GROUP=1
LLM_GROUP1_PROVIDER=anthropic
LLM_GROUP1_BASE_URL=http://127.0.0.1:8082
LLM_GROUP1_API_KEY=key1
LLM_GROUP1_MODEL_HAIKU=claude-haiku-4-5
LLM_GROUP2_PROVIDER=openai
LLM_GROUP2_BASE_URL=https://api.openai.com
LLM_GROUP2_API_KEY=key2
LLM_GROUP2_MODEL_HAIKU=gpt-4o-mini
"""
        root = self._make_env(tmp_path, env_content)
        api_keys, base_urls, model_list = build_mllmcelltype_config(project_root=root)

        assert api_keys == {"anthropic": "key1", "openai": "key2"}
        assert len(base_urls) == 2
        # 无 "#" 后缀
        assert all("#" not in k for k in api_keys)
        assert all("#" not in k for k in base_urls)
        assert set(model_list) == {"claude-haiku-4-5", "gpt-4o-mini"}

    # ------------------------------------------------------------------
    # 守卫：空 key group 不占用路由键
    # ------------------------------------------------------------------

    def test_empty_key_group_does_not_steal_bare_name(self, tmp_path):
        """group1 anthropic 但 key 为空 + 配了 haiku → 不进 api_keys，
        group2 anthropic 有 key → 落在裸名 'anthropic' 上，不被挤到 #2。
        """
        env_content = """LLM_DEFAULT_GROUP=1
LLM_GROUP1_PROVIDER=anthropic
LLM_GROUP1_BASE_URL=http://127.0.0.1:8082
LLM_GROUP1_API_KEY=
LLM_GROUP1_MODEL_HAIKU=claude-haiku-4-5
LLM_GROUP2_PROVIDER=anthropic
LLM_GROUP2_BASE_URL=http://127.0.0.1:8083
LLM_GROUP2_API_KEY=key2
LLM_GROUP2_MODEL_HAIKU=claude-haiku-4-5
"""
        root = self._make_env(tmp_path, env_content)
        api_keys, base_urls, model_list = build_mllmcelltype_config(project_root=root)

        # 只有 group2 的 key 进入（group1 空 key 不贡献）
        assert api_keys == {"anthropic": "key2"}
        assert len(api_keys) == 1
        assert "anthropic" in base_urls
