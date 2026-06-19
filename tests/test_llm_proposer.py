"""Tests for notebooks/_llm_proposer.py — obs 对齐 LLM 提议器纯函数。

本测试文件测四个确定性纯函数，不调 LLM（LLM 调用在 notebook 由 PI 手动跑）。
测试范围：
- build_proposal_prompt: prompt 含关键信息（obs 列名、manifest、本体引用）
- parse_proposal: 正常 JSON / 带代码围栏 / 缺字段容错 / 空输入
- merge_into_manifest: 合并不丢原字段 / 不覆盖已存在映射 / 按字段合并
- write_manifest: round-trip 读回一致 + 中文保留
"""

import importlib.util
import json
import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest
import yaml


# ---------------------------------------------------------------------------
# 加载 notebooks/_llm_proposer（不在 src/ 下，用 importlib 从路径加载）
# ---------------------------------------------------------------------------

def _load_module():
    """从绝对路径加载 notebooks/_llm_proposer.py。"""
    project_root = Path(__file__).resolve().parent.parent
    module_path = project_root / "notebooks" / "_llm_proposer.py"
    spec = importlib.util.spec_from_file_location(
        "_llm_proposer", str(module_path)
    )
    assert spec is not None, f"模块路径不存在: {module_path}"
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_llm = _load_module()


# ---------------------------------------------------------------------------
# 共享 fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_obs_head():
    """模拟 Kim 数据集的 obs.head()"""
    return pd.DataFrame({
        "condition": ["CN", "CN", "Com", "Incom", "Incom",
                       "na", "na", "SI", "Com", "CN"],
        "source_file": ["kim_sample_1"] * 10,
        "_batch": ["batch_A"] * 5 + ["batch_B"] * 5,
        "n_genes": [2500, 3100, 2800, 2200, 3400,
                     2900, 2600, 3100, 2400, 2700],
    })


@pytest.fixture
def sample_manifest():
    """模拟 Kim manifest"""
    return {
        "species": "human",
        "input": {
            "format": "h5ad",
            "path": "data/_subset/kim/kim_subset.h5ad",
            "gene_id_format": "symbol",
        },
        "source_dataset": "Kim_2023",
        "project_id": "gcpl_gastric_2026",
        "disease_system": "gastric",
        "obs_mapping": {
            "sample_id": "source_file",
            "batch": "_batch",
        },
        "original_annotations": [],
    }


@pytest.fixture
def sample_ontology():
    """模拟 gastric.yaml 的 ontology dict"""
    return {
        "ontology": "gastric",
        "nodes": [
            {"id": "normal", "label": "Normal gastric mucosa",
             "parent": None, "mondo": "MONDO:0024516"},
            {"id": "gastric_cancer", "label": "Gastric adenocarcinoma",
             "parent": "dysplasia", "mondo": "MONDO:0001056"},
            {"id": "CAG", "label": "Chronic atrophic gastritis",
             "parent": "chronic_gastritis", "mondo": "MONDO:0005048"},
        ],
    }


@pytest.fixture
def sample_clinical_head():
    """模拟临床信息表头部"""
    return pd.DataFrame({
        "patient_id": ["P001", "P002", "P003"],
        "age": [58, 62, 45],
        "sex": ["M", "F", "M"],
        "H_pylori": ["positive", "negative", "positive"],
    })


# ---------------------------------------------------------------------------
# build_proposal_prompt
# ---------------------------------------------------------------------------


class TestBuildProposalPrompt:
    def test_returns_two_strings(self, sample_obs_head, sample_manifest):
        system, user = _llm.build_proposal_prompt(sample_obs_head, sample_manifest)
        assert isinstance(system, str) and len(system) > 0
        assert isinstance(user, str) and len(user) > 0

    def test_user_contains_obs_column_names(self, sample_obs_head, sample_manifest):
        _, user = _llm.build_proposal_prompt(sample_obs_head, sample_manifest)
        assert "condition" in user
        assert "source_file" in user
        assert "_batch" in user

    def test_user_contains_manifest_info(self, sample_obs_head, sample_manifest):
        _, user = _llm.build_proposal_prompt(sample_obs_head, sample_manifest)
        assert "Kim_2023" in user
        assert "gcpl_gastric_2026" in user

    def test_user_contains_obs_values(self, sample_obs_head, sample_manifest):
        _, user = _llm.build_proposal_prompt(sample_obs_head, sample_manifest)
        # obs 的实际取值出现在 prompt 中
        assert "CN" in user
        assert "Com" in user

    def test_user_contains_ontology_nodes(self, sample_obs_head, sample_manifest,
                                          sample_ontology):
        _, user = _llm.build_proposal_prompt(
            sample_obs_head, sample_manifest,
            ontology_dict=sample_ontology,
        )
        assert "MONDO:0024516" in user
        assert "Normal gastric mucosa" in user
        assert "gastric_cancer" in user

    def test_user_mentions_clinical_table(self, sample_obs_head, sample_manifest,
                                          sample_clinical_head):
        _, user = _llm.build_proposal_prompt(
            sample_obs_head, sample_manifest,
            clinical_head_df=sample_clinical_head,
        )
        assert "patient_id" in user
        assert "H_pylori" in user
        assert "P001" in user

    def test_output_schema_in_user_prompt(self, sample_obs_head, sample_manifest):
        _, user = _llm.build_proposal_prompt(sample_obs_head, sample_manifest)
        assert "obs_mapping" in user
        assert "value_mapping" in user
        assert "ontology" in user
        assert "rationale" in user

    def test_empty_manifest_handled(self, sample_obs_head):
        system, user = _llm.build_proposal_prompt(sample_obs_head, {})
        assert "（manifest 为空" in user


# ---------------------------------------------------------------------------
# parse_proposal
# ---------------------------------------------------------------------------


class TestParseProposal:
    def test_parses_valid_json(self):
        proposal = {
            "obs_mapping": {"disease": "condition", "sex": "Sex"},
            "value_mapping": {"disease": {"CN": "normal"}},
            "ontology": {"disease_ontology_term_id": "MONDO:0005048"},
            "rationale": "基于列名语义匹配。",
        }
        raw = json.dumps(proposal, ensure_ascii=False)
        result = _llm.parse_proposal(raw)
        assert result["obs_mapping"] == {"disease": "condition", "sex": "Sex"}
        assert result["value_mapping"] == {"disease": {"CN": "normal"}}
        assert result["ontology"] == {"disease_ontology_term_id": "MONDO:0005048"}
        assert "基于列名" in result["rationale"]

    def test_strips_markdown_fence(self):
        proposal = {
            "obs_mapping": {"disease": "condition"},
            "value_mapping": {},
            "ontology": {},
            "rationale": "test",
        }
        raw = "```json\n" + json.dumps(proposal, ensure_ascii=False) + "\n```"
        result = _llm.parse_proposal(raw)
        assert result["obs_mapping"] == {"disease": "condition"}

    def test_strips_markdown_fence_no_json_tag(self):
        proposal = {
            "obs_mapping": {"disease": "condition"},
            "value_mapping": {},
            "ontology": {},
            "rationale": "test",
        }
        raw = "```\n" + json.dumps(proposal, ensure_ascii=False) + "\n```"
        result = _llm.parse_proposal(raw)
        assert result["obs_mapping"] == {"disease": "condition"}

    def test_missing_fields_fallback_to_empty(self):
        raw = json.dumps({"obs_mapping": {"disease": "condition"}})
        result = _llm.parse_proposal(raw)
        assert result["obs_mapping"] == {"disease": "condition"}
        assert result["value_mapping"] == {}
        assert result["ontology"] == {}
        assert result["rationale"] == ""

    def test_empty_string_returns_fallback(self):
        result = _llm.parse_proposal("")
        assert result["obs_mapping"] == {}
        assert result["value_mapping"] == {}
        assert result["ontology"] == {}
        assert result["rationale"] == ""

    def test_none_text_returns_fallback(self):
        result = _llm.parse_proposal(None)
        assert result == {
            "obs_mapping": {},
            "value_mapping": {},
            "ontology": {},
            "rationale": "",
        }

    def test_invalid_json_returns_fallback(self):
        result = _llm.parse_proposal("not even json at all")
        assert result["obs_mapping"] == {}

    def test_trailing_comma_recovery(self):
        # LLM 常见错误：对象末尾多一个逗号
        raw = '{"obs_mapping": {"disease": "condition",}, "value_mapping": {}, "ontology": {}, "rationale": "ok"}'
        result = _llm.parse_proposal(raw)
        assert result["obs_mapping"] == {"disease": "condition"}

    def test_non_dict_values_coerced(self):
        raw = json.dumps({
            "obs_mapping": ["not", "a", "dict"],  # 错误类型
            "value_mapping": {},
            "ontology": {},
            "rationale": 12345,  # 非字符串
        })
        result = _llm.parse_proposal(raw)
        assert result["obs_mapping"] == {}
        assert result["rationale"] == "12345"

    def test_extracts_json_from_mixed_text(self):
        # LLM 有时在 JSON 前后加解释文字
        proposal = {
            "obs_mapping": {"disease": "condition"},
            "value_mapping": {},
            "ontology": {},
            "rationale": "test",
        }
        raw = "这是分析结果：\n" + json.dumps(proposal, ensure_ascii=False) + "\n以上就是我的建议。"
        result = _llm.parse_proposal(raw)
        assert result["obs_mapping"] == {"disease": "condition"}


# ---------------------------------------------------------------------------
# merge_into_manifest
# ---------------------------------------------------------------------------


class TestMergeIntoManifest:
    def test_adds_new_obs_mapping(self, sample_manifest):
        proposal = {
            "obs_mapping": {"disease": "condition"},
            "value_mapping": {},
            "ontology": {},
            "rationale": "",
        }
        result = _llm.merge_into_manifest(sample_manifest, proposal)
        assert result["obs_mapping"]["disease"] == "condition"
        # 原有映射保留
        assert result["obs_mapping"]["sample_id"] == "source_file"
        assert result["obs_mapping"]["batch"] == "_batch"

    def test_does_not_overwrite_existing_obs_mapping(self, sample_manifest):
        proposal = {
            "obs_mapping": {
                "sample_id": "different_column",  # 尝试覆盖
                "disease": "condition",           # 新增
            },
            "value_mapping": {},
            "ontology": {},
            "rationale": "",
        }
        result = _llm.merge_into_manifest(sample_manifest, proposal)
        # 已有映射不被覆盖
        assert result["obs_mapping"]["sample_id"] == "source_file"
        # 新映射正常添加
        assert result["obs_mapping"]["disease"] == "condition"

    def test_merges_value_mapping_incrementally(self, sample_manifest):
        # 先给 manifest 加一个已有的 value_mapping
        manifest = dict(sample_manifest)
        manifest["value_mapping"] = {
            "disease": {"CN": "normal"},
        }
        proposal = {
            "obs_mapping": {},
            "value_mapping": {
                "disease": {
                    "CN": "WRONG",      # 尝试覆盖已有键
                    "Com": "complete_IM",  # 新键
                },
                "sex": {"M": "male"},  # 新字段
            },
            "ontology": {},
            "rationale": "",
        }
        result = _llm.merge_into_manifest(manifest, proposal)
        assert result["value_mapping"]["disease"]["CN"] == "normal"  # 不覆盖
        assert result["value_mapping"]["disease"]["Com"] == "complete_IM"  # 新增
        assert result["value_mapping"]["sex"]["M"] == "male"  # 新字段

    def test_fills_empty_ontology_fields(self, sample_manifest):
        manifest = dict(sample_manifest)
        manifest["ontology"] = {"disease_ontology_term_id": ""}
        proposal = {
            "obs_mapping": {},
            "value_mapping": {},
            "ontology": {
                "disease_ontology_term_id": "MONDO:0005048",
                "tissue_ontology_term_id": "UBERON:0001199",
            },
            "rationale": "",
        }
        result = _llm.merge_into_manifest(manifest, proposal)
        assert result["ontology"]["disease_ontology_term_id"] == "MONDO:0005048"
        assert result["ontology"]["tissue_ontology_term_id"] == "UBERON:0001199"

    def test_does_not_overwrite_existing_ontology(self, sample_manifest):
        manifest = dict(sample_manifest)
        manifest["ontology"] = {"disease_ontology_term_id": "MONDO:0001056"}
        proposal = {
            "obs_mapping": {},
            "value_mapping": {},
            "ontology": {
                "disease_ontology_term_id": "MONDO:9999999",  # 尝试覆盖
            },
            "rationale": "",
        }
        result = _llm.merge_into_manifest(manifest, proposal)
        assert result["ontology"]["disease_ontology_term_id"] == "MONDO:0001056"

    def test_preserves_unrelated_fields(self, sample_manifest):
        proposal = {
            "obs_mapping": {"disease": "condition"},
            "value_mapping": {},
            "ontology": {},
            "rationale": "",
        }
        result = _llm.merge_into_manifest(sample_manifest, proposal)
        assert result["species"] == "human"
        assert result["project_id"] == "gcpl_gastric_2026"
        assert result["disease_system"] == "gastric"
        assert result["source_dataset"] == "Kim_2023"
        assert result["input"]["format"] == "h5ad"
        assert result["original_annotations"] == []

    def test_does_not_modify_input(self, sample_manifest):
        original = dict(sample_manifest)
        proposal = {
            "obs_mapping": {"disease": "condition"},
            "value_mapping": {},
            "ontology": {},
            "rationale": "",
        }
        _llm.merge_into_manifest(sample_manifest, proposal)
        # 原始 dict 不变
        assert sample_manifest == original

    def test_non_dict_value_mapping_skipped(self, sample_manifest):
        proposal = {
            "obs_mapping": {},
            "value_mapping": {"disease": "not_a_dict"},
            "ontology": {},
            "rationale": "",
        }
        result = _llm.merge_into_manifest(sample_manifest, proposal)
        # 不应该崩溃，且不写入错误类型的 value
        assert "disease" not in result.get("value_mapping", {})


# ---------------------------------------------------------------------------
# P1 xfail: merge_into_manifest — 已有 value_mapping 条目非 dict 时崩溃
# ---------------------------------------------------------------------------


class TestMergeIntoManifestBug:
    """确认 bug：manifest 中 value_mapping 某条目是 string 而非 dict 时崩溃。"""

    @pytest.mark.xfail(
        reason=(
            "BUG: merge_into_manifest 未检查 existing_val[norm_name] 是否为 dict"
            "（_llm_proposer.py:318）。当手动编辑 YAML 把 value_mapping 某 key"
            "误写为 plain string 而非 dict 时，line 318 'src_val not in "
            "existing_val[norm_name]' 对 string 做子串匹配，line 319"
            " existing_val[norm_name][src_val]=norm_val 导致 TypeError 崩溃。"
            "修复方向：line 315 后加 isinstance(existing_val.get(norm_name), dict) 守卫"
        ),
        strict=True,
    )
    def test_value_mapping_existing_str_not_dict_crashes(self, sample_manifest):
        """manifest 已有 value_mapping.disease='not_a_dict'（string）→ 应优雅处理而非崩溃。"""
        manifest = dict(sample_manifest)
        manifest["value_mapping"] = {"disease": "not_a_dict"}  # 手动 YAML 编辑误写
        proposal = {
            "obs_mapping": {},
            "value_mapping": {
                "disease": {"CAG": "atrophic_gastritis"},
            },
            "ontology": {},
            "rationale": "",
        }
        # 期望：不崩溃，优雅跳过非 dict 条目
        # 当前行为：TypeError 崩溃
        result = _llm.merge_into_manifest(manifest, proposal)
        # 即使 manifest 原值为 string，也不应崩溃
        assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# write_manifest + round-trip
# ---------------------------------------------------------------------------


class TestWriteManifest:
    def test_round_trip_preserves_data(self, sample_manifest):
        """写回后 yaml.safe_load 读回，关键字段一致。"""
        with tempfile.NamedTemporaryFile(
            suffix=".yaml", mode="w", delete=False, encoding="utf-8"
        ) as f:
            tmp_path = f.name

        try:
            _llm.write_manifest(sample_manifest, tmp_path)
            with open(tmp_path, "r", encoding="utf-8") as f:
                reloaded = yaml.safe_load(f)
            assert reloaded["species"] == "human"
            assert reloaded["source_dataset"] == "Kim_2023"
            assert reloaded["disease_system"] == "gastric"
            assert reloaded["obs_mapping"]["sample_id"] == "source_file"
            assert reloaded["obs_mapping"]["batch"] == "_batch"
        finally:
            os.unlink(tmp_path)

    def test_chinese_preserved(self):
        """中文内容写回后不变成 Unicode 转义序列。"""
        manifest = {"description": "慢性萎缩性胃炎患者胃粘膜活检标本"}
        with tempfile.NamedTemporaryFile(
            suffix=".yaml", mode="w", delete=False, encoding="utf-8"
        ) as f:
            tmp_path = f.name

        try:
            _llm.write_manifest(manifest, tmp_path)
            raw_text = Path(tmp_path).read_text(encoding="utf-8")
            # 中文字符原样存在，不是 \\uXXXX 转义
            assert "慢性萎缩性胃炎" in raw_text
            assert "\\u" not in raw_text
        finally:
            os.unlink(tmp_path)

    def test_keys_not_sorted(self):
        """sort_keys=False 确保键按写入顺序排列（或自然顺序），不是字母序。"""
        manifest = {
            "z_field": "last_in_alphabet",
            "a_field": "first_in_alphabet",
            "m_field": "middle",
        }
        with tempfile.NamedTemporaryFile(
            suffix=".yaml", mode="w", delete=False, encoding="utf-8"
        ) as f:
            tmp_path = f.name

        try:
            _llm.write_manifest(manifest, tmp_path)
            raw_text = Path(tmp_path).read_text(encoding="utf-8")
            # z_field 在第一行 "z_field:" 应出现在 a_field 之前
            z_pos = raw_text.find("z_field:")
            a_pos = raw_text.find("a_field:")
            assert z_pos < a_pos, (
                f"键顺序被排序！z_field 在 {z_pos}, a_field 在 {a_pos}"
            )
        finally:
            os.unlink(tmp_path)

    def test_file_ends_with_newline(self):
        """POSIX 惯例：文件末尾有换行。"""
        manifest = {"key": "value"}
        with tempfile.NamedTemporaryFile(
            suffix=".yaml", mode="w", delete=False, encoding="utf-8"
        ) as f:
            tmp_path = f.name

        try:
            _llm.write_manifest(manifest, tmp_path)
            raw_text = Path(tmp_path).read_text(encoding="utf-8")
            assert raw_text.endswith("\n")
        finally:
            os.unlink(tmp_path)
