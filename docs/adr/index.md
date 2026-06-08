---
title: "ADR 索引"
type: adr-index
project_id: "scrna-integration-framework"
created: "2026-06-08"
updated: "2026-06-08"
---

# 架构决策记录（ADR）索引

本目录记录 scRNA-seq 整合分析框架在设计与实现过程中做出的关键架构决策。
每份 ADR 包含问题背景、候选方案、决策理由与后果。

ADR 编号按决策时间排序，非按重要程度。

| 编号 | 标题 | 状态 | 一句话摘要 |
|------|------|------|-----------|
| [ADR-0001](0001-thin-framework-over-scanpy.md) | Thin framework over native scanpy | accepted | 框架不包装 scanpy 函数，直接使用 scanpy 原生 API；仅在 scanpy 覆盖不到的 4 个缺口（多源 IO、参数扫描、阶段报告、运行追踪）添加薄层 |
| [ADR-0002](0002-r-bridge-rpy2.md) | R interoperability via rpy2 + anndata2ri | **superseded** by ADR-0007 | 最初计划通过 rpy2 统一处理所有 R 互操作；后继因重型 R 工具的内存与版本冲突问题，被 ADR-0007 的"按工具分流"策略取代 |
| [ADR-0003](0003-plain-code-over-plugin-systems.md) | Plain code over plugin systems | accepted | 框架接受用户扩展时，默认方案是普通函数（或函数列表）作为参数传入；拒绝装饰器注册、插件发现、动态导入等"框架魔法" |
| [ADR-0004](0004-framework-deletion-log.md) | Framework deletion log: from 11 si.* APIs to 2 functions | accepted | PI 第二次审查将框架从 11 个公开 API 砍到 2 个函数（后经 ADR-0005 增至 3 个），拒绝一切包装 scanpy/anndata/pandas/yaml 既有功能的抽象层 |
| [ADR-0005](0005-load-markers-helper.md) | `load_markers` as a legitimate third framework function | accepted | 将标记物加载 + role 过滤 + 分组操作集中为一个函数，作为 ADR-0004 逃生舱门槛的首次触发案例；防止 `negative` 标记物被误用于基因集评分 |
| [ADR-0006](0006-yaml-manifest-over-notebook-config.md) | YAML manifest for dataset facts | accepted | 每个数据集的"事实"（来源路径、批次、物种、已知细胞类型）放在 `manifest.yaml` 中，而非 notebook PARAMS cell 的 Python dict；YAML 只存不变事实，不存参数 |
| [ADR-0007](0007-r-bridge-tool-split.md) | R bridge split by tool | accepted (supersedes ADR-0002) | R 互操作按工具分三路：纯 Python 替代优先、rpy2 用于轻量 notebook 内调用、subprocess `Rscript` 用于重型 R 工具（Monocle3、hdWGCNA） |
| [ADR-0008](0008-absorb-student-code-by-rewriting.md) | Absorb student-code downstream techniques by re-implementing | accepted | student-code 中的下游技术（entropy、CytoTRACE、UCell、scCODA 等）全部吸收，但每个都按本项目规范重新实现，不直接复制——student code 是"算法正确性参考"而非代码捐献者 |
| [ADR-0009](0009-de-encapsulation-teaching-transparency.md) | 教学透明：删 sweep、摊平 IO 琐碎 helper、注释中文化 | accepted | 从薄框架进一步回退到"教学透明"：删除 `sweep` 函数、摊平 IO 中的琐碎 helper、所有注释改为中文；判据是"非 CS 学生打开 notebook 能否逐行看懂" |

## 引用 ADR

- 代码注释或 PR 描述中引用 ADR 时使用格式：`ADR-XXXX`（如 `ADR-0005`）
- 新 ADR 提案应在 PR 描述中明确引用将被修改/替代的已有 ADR
- Superseded 的 ADR 不删除，保留在目录中并在状态栏标注被哪个 ADR 替代

## 相关文档

- [项目 SPEC](../SPEC.md) — 实现规格，反映当前 ADR 的落地状态
- [项目 CONTEXT](../CONTEXT.md) — 术语表与架构姿态
- [框架根 README](../../README.md)
