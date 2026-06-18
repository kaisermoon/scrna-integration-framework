---
status: accepted
---

# ADR-0014: obs 对齐两相设计——LLM 设计期提议、运行期纯确定性，向 CellxGene 字段标准对齐

- **Status**: Accepted
- **Date**: 2026-06-18
- **Supersedes**: None
- **Related**: ADR-0006 (YAML manifest 作为数据事实单一源), ADR-0009 (去封装与教学透明), ADR-0004 (薄框架), SPEC「obs Schema」节

## Context

跨数据集 obs 列高度异构：同一语义（如疾病状态）在各源用不同列名（Kim `condition` / Nowicki `Patient_status` / Nancang `disease_state`）、不同取值词表（`GC` / `gastric cancer` / `Incom`），部分源是隐晦编码（Kim 的 `Com`/`Incom`/`na`），CellxGene 七字段（`disease` / `disease_ontology_term_id` / `tissue` / `tissue_ontology_term_id` / `assay` / `sex` / `development_stage`）常缺失或为自由文本。

2026-06-18 对 `read_with_manifest` 的 12 场景压力测试确认对齐机制当前**不能可靠对齐 Core obs schema**：6 个 P0 代码缺陷（Layer 1 完整性零校验、`project_id` 从未写入 obs、obs_mapping 目标列静默覆盖、NaN 经 `astype(str)` 不可逆转为 `"nan"`、value_mapping 未覆盖值静默穿透、project_specific 源列缺失静默跳过）+ 4 个 P1 设计缺口。问题全部出在**确定性应用层缺校验/缺警告**，无一与 LLM 相关。

PI 提出"真正的对齐是否必须依赖 LLM"。厘清后达成两点定调：① 认同两相设计；② obs 字段尽可能向 CellxGene 标准对齐；③ 本体接地用 LLM 提议。

## Decision

**把"对齐"拆成两个根本不同的动作，用一条冻结边界隔开：**

```
设计期（LLM 提议 + PI 确认）  →  冻结为 manifest.yaml  →  运行期（纯确定性，零 LLM）
```

### 一、运行期：纯确定性，永不调 LLM

`read_with_manifest` 是确定性函数。manifest 是数据事实单一源（ADR-0006），运行时只**应用**已冻结的映射，严格校验，绝不调 LLM。理由：可复现、可审计（映射可 diff）、可教学透明（ADR-0009）、772k 细胞规模成本可控、科学结论的判断权不外包（SOUL.md 底线）。LLM 提议的映射是**假设非真相**，真相来自数据集原文 + PI 领域知识，故 LLM 产物必须经 PI 确认并冻结后才进运行期。

运行期对齐的正确行为（修复 P0/P1 的判据）：

1. **Layer 1 完整性 fail loudly**：`source_dataset` / `project_id` / `disease_system` 缺失即抛错；三者均写入 obs（`project_id` 当前缺失，必补）。
2. **obs_mapping 不静默**：源列缺失 warn、目标列已存在冲突 warn（不静默覆盖）、NaN 保留为 NaN（不转 `"nan"` 字符串）。
3. **value_mapping 覆盖率可见**：未覆盖值不静默穿透，发警告列出未映射取值。
4. **行为一致**：project_specific 源列缺失补 warn，与 obs_mapping / original_annotations 对齐。

### 二、设计期：LLM 提议器，PI 确认，写回 manifest

新数据集首次接入时，一个**设计期工具**（notebook cell 或独立 helper，非 `read_with_manifest` 运行路径）读 obs 头部 + manifest + 临床表头部，**提议** obs_mapping / value_mapping / CellxGene 七字段值 / 本体接地（自由文本 → MONDO / UBERON term ID）。PI 在 notebook 内逐条确认，确认后写回 `manifest.yaml`。这正是 SPEC「Layer 2」节已描述、因未配 key 而搁置的 best-effort fix，本 ADR 将其正式定为架构原则。

LLM 必要性分级（不均匀，避免过度依赖）：

| 子问题 | LLM 角色 | 权威来源 |
|---|---|---|
| 清晰列名 → 规范字段、清晰取值统一 | 加速，非必需（人眼可判） | PI 看一次 obs |
| 隐晦编码解码、CellxGene 字段填充 | 提议（有幻觉风险） | 数据集原文 + PI |
| 本体接地（文本 → MONDO/UBERON ID） | **提议（PI 选定 LLM 方案）** | LLM 提议 + PI 确认 |
| 未映射值检测 / 列冲突检测 | **不参与** | 确定性代码 |

### 三、向 CellxGene 字段标准对齐

设计期提议优先把 obs 对齐到 CellxGene schema 的七个 Layer 2 字段及其 ontology term ID。这使数据集对外可互操作、对内可跨源整合（DEG / abundance / 注释 / 跨组织比较均依赖这些字段）。Layer 3 项目自定义字段（`H_pylori_status` 等）原样透传。

## Consequences

**正向**：对齐可规模化（LLM 减轻 PI 逐源手写 manifest 的负担）同时保持可信（人确认 + 冻结）；运行期可复现可审计；P0/P1 修复有明确判据；与既有 ADR（薄框架、去封装、manifest 单一源）一致。

**代价 / 边界**：设计期 LLM 需 OpenRouter key（PI 待办，与 stage6 共用）；本体接地可能引入 LLM 调用成本；LLM 提议质量依赖 obs 头部信息量，信息不足的字段保持 NaN 并 warn，不强行猜测。

**不做**：全自动 LLM 运行时对齐（不可复现 / 不可审计 / 幻觉污染科学结论 / 成本不可控 / 判断权外包）；运行期任何 LLM 调用。

## 落地顺序

1. **先做（与 LLM 无关，地基）**：修复 6 个 P0 + 4 个 P1，落实运行期确定性对齐四条行为。详见 `results/stress_test_obs_alignment_2026-06-18.md`。
2. **后做（需 key）**：设计期 LLM 提议器实现（SPEC Layer 2 best-effort fix），含本体接地。

> 同步：SPEC「obs Schema」节补"两相边界"说明；`_plan.md` 新增对齐修复 PR。
