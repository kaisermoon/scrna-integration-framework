"""设计期 obs 对齐 LLM 提议器——ADR-0014 第二相。

本模块是**确定性纯函数**集合，不调 LLM。LLM 调用在 notebook 里完成。
每个函数可独立单测，不依赖 src/scrna_integration/ 框架（只读不写）。

四个公开函数：
- ``build_proposal_prompt()`` — 构造给 LLM 的结构化 prompt
- ``parse_proposal()`` — 解析 LLM 返回的 JSON 提议（容错）
- ``merge_into_manifest()`` — 把 PI 确认后的提议合并进 manifest
- ``write_manifest()`` — YAML 写回，保留中文与顺序

面向非计算机专业 PI/学生：从上到下线性可读。
"""

from __future__ import annotations

import copy
import json
import re
from typing import Any

import yaml


# ---------------------------------------------------------------------------
# 1. build_proposal_prompt
# ---------------------------------------------------------------------------


def build_proposal_prompt(
    obs_head_df: "pandas.DataFrame",
    manifest_dict: dict[str, Any],
    clinical_head_df: "pandas.DataFrame | None" = None,
    ontology_dict: dict[str, Any] | None = None,
    disease_system: str = "gastric",
) -> tuple[str, str]:
    """构造给 LLM 的 obs 对齐提议 prompt。

    把 obs 列名 + dtype + 前几行取值、现有 manifest、临床表头、本体参考
    塞进 prompt，要求 LLM 返回结构化 JSON。

    Parameters
    ----------
    obs_head_df : pd.DataFrame
        adata.obs.head(10) 的输出，包含列名、dtype、前几行实际值。
    manifest_dict : dict
        现有 manifest YAML 的 Python dict（可为空或部分填充）。
    clinical_head_df : pd.DataFrame or None
        临床信息表头部（如有）。None 表示无临床表。
    ontology_dict : dict or None
        疾病本体参考（如 gastric.yaml 的 dict）。None 表示无本体参考。
    disease_system : str
        疾病系统标识，如 "gastric"、"synovium"。

    Returns
    -------
    (system, user) : tuple[str, str]
        system prompt 和 user prompt，可直接送给 LLM。
    """
    # ----- system prompt：角色 + 输出格式约束 -----
    system = (
        f"你是单细胞转录组学数据集成专家，专精于跨数据集的 obs 字段对齐。\n"
        f"当前疾病系统：{disease_system}。\n\n"
        f"你的任务是：阅读一个数据集的 obs 列名与取值、现有 manifest、"
        f"临床表头与本体参考，**提议** obs_mapping / value_mapping / ontology 三块。\n\n"
        f"你的提议是**假设而非真相**——PI 会逐条确认后才写入 manifest。\n"
        f"对于无法判断的字段，明确标注 unknown 或留空，不要强行猜测。"
    )

    # ----- user prompt：结构化信息注入 -----
    lines: list[str] = []

    # 1. obs 列信息
    lines.append("## 1. obs 列信息（adata.obs.head(10)）\n")
    lines.append("| 列名 | dtype | 前 5 个非空值 |")
    lines.append("|------|-------|--------------|")
    for col in obs_head_df.columns:
        dtype = str(obs_head_df[col].dtype)
        vals = obs_head_df[col].dropna().head(5).astype(str).tolist()
        vals_str = ", ".join(vals) if vals else "（全空）"
        lines.append(f"| {col} | {dtype} | {vals_str} |")
    lines.append("")

    # 2. 现有 manifest
    lines.append("## 2. 现有 manifest\n")
    if manifest_dict:
        lines.append("```yaml")
        # 用 safe_dump 展示当前 manifest
        lines.append(yaml.safe_dump(manifest_dict, allow_unicode=True, sort_keys=False).rstrip())
        lines.append("```")
    else:
        lines.append("（manifest 为空——新数据集首次接入）")
    lines.append("")

    # 3. 临床表头部
    lines.append("## 3. 临床信息表头部\n")
    if clinical_head_df is not None and not clinical_head_df.empty:
        lines.append("| 列名 | dtype | 前 3 个非空值 |")
        lines.append("|------|-------|--------------|")
        for col in clinical_head_df.columns:
            dtype = str(clinical_head_df[col].dtype)
            vals = clinical_head_df[col].dropna().head(3).astype(str).tolist()
            vals_str = ", ".join(vals) if vals else "（全空）"
            lines.append(f"| {col} | {dtype} | {vals_str} |")
    else:
        lines.append("（无临床信息表）")
    lines.append("")

    # 4. 本体参考
    lines.append("## 4. 本体参考\n")
    if ontology_dict:
        lines.append("以下为已知疾病本体的 id / label / MONDO 映射：\n")
        nodes = ontology_dict.get("nodes", [])
        if nodes:
            for node in nodes:
                nid = node.get("id", "?")
                label = node.get("label", "?")
                mondo = node.get("mondo", "")
                mondo_str = f" (MONDO:{mondo})" if mondo else ""
                lines.append(f"- `{nid}`: {label}{mondo_str}")
        lines.append("\n请优先使用上述 MONDO term ID 做疾病本体接地。"
                       "如 obs 中出现的疾病在上述本体中找不到精确匹配，"
                       "尝试找最近祖先节点。")
    else:
        lines.append("（无本体参考——LLM 可据领域知识提议 MONDO/UBERON term ID）")
    lines.append("")

    # 5. 输出格式约束
    lines.append("## 5. 请返回结构化 JSON\n")
    lines.append("```json")
    lines.append("{")
    lines.append('  "obs_mapping": {')
    lines.append('    // 源列名 → 规范字段名。规范字段名只能是:')
    lines.append('    //   sample_id / donor_id / batch / disease / disease_ontology_term_id')
    lines.append('    //   / tissue / tissue_ontology_term_id / assay / sex / development_stage')
    lines.append('    // 无法映射的列不写入；无对应列则字段留空 {}')
    lines.append('    "disease": "condition",        // 例：源列 condition → 规范字段 disease')
    lines.append('    "sex": "Sex"                    // 例：源列 Sex → 规范字段 sex')
    lines.append("  },")
    lines.append('  "value_mapping": {')
    lines.append('    // 对需要归一化的字段，列出 源取值 → 规范取值 的映射')
    lines.append('    // 无法判断的取值不写入')
    lines.append('    "disease": {')
    lines.append('      "CN": "normal",              // 隐晦编码 → 可读标签')
    lines.append('      "GC": "gastric_cancer"')
    lines.append("    },")
    lines.append('    "sex": {')
    lines.append('      "M": "male",')
    lines.append('      "F": "female"')
    lines.append("    }")
    lines.append("  },")
    lines.append('  "ontology": {')
    lines.append('    // CellxGene 七字段的本体接地（直接填 term ID 字符串）')
    lines.append('    "disease_ontology_term_id": "MONDO:0005048",  // 慢性萎缩性胃炎')
    lines.append('    "tissue_ontology_term_id": "UBERON:0001199",   // 胃粘膜')
    lines.append('    // 上述两项自动从 value_mapping 或 obs 中推断。')
    lines.append('    // 其他五字段根据 obs 现有信息判断，无法判断留空字符串')
    lines.append('    "assay": "EFO:0008913",         // 或留空 ""')
    lines.append('    "sex": "PATO:0000384",')
    lines.append('    "development_stage": ""')
    lines.append("  },")
    lines.append('  "rationale": "简要解释你的推断依据（≤200 字），标注不确定处"')
    lines.append("}")
    lines.append("```")
    lines.append("\n**IMPORTANT**: 直接返回 JSON，不要用 markdown 代码围栏包裹。")

    user = "\n".join(lines)
    return system, user


# ---------------------------------------------------------------------------
# 2. parse_proposal
# ---------------------------------------------------------------------------


def parse_proposal(llm_text: str) -> dict[str, Any]:
    """解析 LLM 返回的文本，提取 obs 对齐提议 JSON。

    容错处理：
    - 剥离 `` ```json ... ``` `` 代码围栏
    - 缺 ``obs_mapping`` / ``value_mapping`` / ``ontology`` 字段时给空块
    - 保底返回空块 dict（不含 ``rationale`` 则给空字符串）

    Parameters
    ----------
    llm_text : str
        LLM 原始返回文本。

    Returns
    -------
    dict
        ``{obs_mapping, value_mapping, ontology, rationale}``。
        即使解析失败也返回含空块的 dict，不抛异常。
    """
    fallback: dict[str, Any] = {
        "obs_mapping": {},
        "value_mapping": {},
        "ontology": {},
        "rationale": "",
    }

    if not llm_text or not llm_text.strip():
        return fallback

    # 策略 1：匹配 ```json ... ``` 代码块
    candidate: str | None = None
    match = re.search(r"```(?:json)?\s*\n?(.*?)\n?```", llm_text, re.DOTALL)
    if match:
        candidate = match.group(1).strip()
    else:
        # 策略 2：找第一个 { 和最后一个 } 之间的内容
        start = llm_text.find("{")
        end = llm_text.rfind("}")
        if start >= 0 and end > start:
            candidate = llm_text[start:end + 1].strip()
        else:
            return fallback

    # 解析 JSON
    try:
        parsed = json.loads(candidate)
    except (json.JSONDecodeError, ValueError):
        # 策略 3：尝试修复常见错误——尾逗号、单引号等
        try:
            # 简单修复：移除尾逗号
            fixed = re.sub(r",\s*([}\]])", r"\1", candidate)
            parsed = json.loads(fixed)
        except (json.JSONDecodeError, ValueError):
            return fallback

    if not isinstance(parsed, dict):
        return fallback

    # 提取四个字段，缺字段给空块
    obs_mapping = parsed.get("obs_mapping", {})
    value_mapping = parsed.get("value_mapping", {})
    ontology = parsed.get("ontology", {})
    rationale = parsed.get("rationale", "")

    # 类型保底
    if not isinstance(obs_mapping, dict):
        obs_mapping = {}
    if not isinstance(value_mapping, dict):
        value_mapping = {}
    if not isinstance(ontology, dict):
        ontology = {}
    if not isinstance(rationale, str):
        rationale = str(rationale)

    return {
        "obs_mapping": obs_mapping,
        "value_mapping": value_mapping,
        "ontology": ontology,
        "rationale": rationale,
    }


# ---------------------------------------------------------------------------
# 3. merge_into_manifest
# ---------------------------------------------------------------------------


def merge_into_manifest(
    manifest_dict: dict[str, Any],
    confirmed_proposal: dict[str, Any],
) -> dict[str, Any]:
    """把 PI 确认后的提议合并进 manifest，返回新 dict。

    合并规则（ADR-0014 设计期语义）：
    - ``obs_mapping``：对于 proposal 中的每个 key（规范字段名），
      如果 manifest 中已存在同名映射且值与 proposal 不同，
      **不覆盖**——保留 manifest 现有值（只新增不存在的 key）。
    - ``value_mapping``：对于 proposal 中的每个字段（如 ``disease``），
      如果 manifest 中已有该字段的映射，**合并**（只新增不存在的取值 key；
      已存在的取值不覆盖）。
    - ``ontology``：如果 manifest 中尚无 ``ontology`` 块，
      直接写入。如果已有，只填充 proposal 中有值而 manifest 中为空的 key。
    - ``rationale``：不写入 manifest（仅供 PI 查阅，存 uns 或丢弃）。
    - 保留 manifest 中所有其他字段（``species``、``input``、
      ``source_dataset``、``project_id``、``disease_system``、
      ``original_annotations``、``preprocessing_done``、``qc_overrides`` 等）。

    Parameters
    ----------
    manifest_dict : dict
        当前 manifest dict（可能为空或部分填充）。
    confirmed_proposal : dict
        PI 确认后的 proposal dict，格式同 ``parse_proposal`` 返回值。

    Returns
    -------
    dict
        合并后的新 manifest dict（深拷贝，不修改输入）。
    """
    result = copy.deepcopy(manifest_dict)

    obs_mapping_new = confirmed_proposal.get("obs_mapping", {})
    value_mapping_new = confirmed_proposal.get("value_mapping", {})
    ontology_new = confirmed_proposal.get("ontology", {})

    # --- obs_mapping: 只新增不存在的 key ---
    if obs_mapping_new:
        existing_obs = result.setdefault("obs_mapping", {})
        for norm_name, src_col in obs_mapping_new.items():
            if norm_name not in existing_obs:
                existing_obs[norm_name] = src_col

    # --- value_mapping: 按字段合并 ---
    if value_mapping_new:
        existing_val = result.setdefault("value_mapping", {})
        for norm_name, val_dict in value_mapping_new.items():
            if not isinstance(val_dict, dict):
                continue
            if norm_name not in existing_val:
                existing_val[norm_name] = {}
            # 合并：只新增不存在的取值
            for src_val, norm_val in val_dict.items():
                if src_val not in existing_val[norm_name]:
                    existing_val[norm_name][src_val] = norm_val

    # --- ontology: 填充空白 ---
    if ontology_new:
        existing_ont = result.setdefault("ontology", {})
        for key, val in ontology_new.items():
            if val and (key not in existing_ont or not existing_ont[key]):
                existing_ont[key] = val

    return result


# ---------------------------------------------------------------------------
# 4. write_manifest
# ---------------------------------------------------------------------------


def write_manifest(manifest_dict: dict[str, Any], path: str) -> None:
    """把 manifest dict 写回 YAML 文件。

    使用 ``yaml.safe_dump`` 确保：
    - 中文原样保留（``allow_unicode=True``）
    - 键不自动排序（``sort_keys=False``）
    - 可读缩进（``indent=2``，block 风格）

    写入后末尾加一个换行（POSIX 惯例）。

    Parameters
    ----------
    manifest_dict : dict
        待写入的 manifest dict。
    path : str
        目标 YAML 文件绝对路径。
    """
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(
            manifest_dict,
            f,
            allow_unicode=True,
            sort_keys=False,
            indent=2,
            default_flow_style=False,
        )
    # POSIX 惯例：文件末尾加换行（safe_dump 默认不带）
    with open(path, "a", encoding="utf-8") as f:
        f.write("\n")
