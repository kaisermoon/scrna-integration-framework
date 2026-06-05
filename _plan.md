---
title: "执行计划：scRNA-seq整合分析框架"
updated: "2026-06-05"
phase: planning
---

## 当前方向

把项目构思打磨为可实施方案，并优先打通"数据接入 + 标准化 + QC"最小可用闭环；首期以 GCPL 现有数据集作端到端验证。

## 任务表

<!-- agent 执行前读此表，执行后更新状态 -->

| # | 任务 | 状态 | 委派 | 输入/约束 | 产出 |
|---|------|------|------|-----------|------|
| 1 | 与 PI 进行 sci-interview，把项目构思蒸馏为可实施方案 | pending | 主 Agent | `planning/项目构思-原始版.md` | `planning/2026-06-05-方案打磨.md` + 更新 `_plan.md` |
| 2 | 决策：项目工程形态（CLI+YAML / Python API / notebook 主导） | pending | 主 Agent | sci-interview 结论 | `_memory.md` 关键决策 |
| 3 | 决策：测试与基准策略（pytest 单元 / 端到端 / scIB benchmark） | pending | 主 Agent | sci-interview 结论 | `_memory.md` 关键决策 |
| 4 | 设计 IO 子系统：cellranger/h5ad/RData reader + CellxGene schema | blocked | coder | 任务 1-3 完成 | `src/io/` 包 + 单元测试 |
| 5 | 设计 QC 子系统：可配置过滤策略 + QC 报告（前后对比图表） | blocked | coder | 任务 4 完成 | `src/qc/` 包 + notebook 02 |
| 6 | 设计降维去批次子系统：PCA/Harmony/scVI/scANVI 多方法对比 | blocked | coder | 任务 5 完成 | `src/embedding/` 包 + notebook 04 |
| 7 | 端到端验证：GCPL 数据跑通 1-6 流水线 | blocked | coder | 任务 6 完成 | `notebooks/end-to-end-pilot.ipynb` |

**状态值**：`pending` | `in_progress` | `done` | `failed` | `blocked`

## 阻塞

- 任务 4-7 全部依赖任务 1-3 的方案打磨，需要 PI 在下一轮对话中把构思讨论清楚
- 大数据存放策略：本地（mac mini / 服务器）/ 学校共享存储 / 云？影响 `data/` 是否做软链 vs 实际拷贝

## 关键决策

| 日期 | 决策 | 理由 |
|------|------|------|
| 2026-06-05 | 项目工程形态采用 `src/{module}/` Python 包 + `notebooks/` 流水线 + `tests/` 单元测试 | 既要可发布、也要可复现；notebook 留住 PI/学生交互打磨能力 |
| 2026-06-05 | 框架以 anndata/scanpy 为内存数据模型主干，CellxGene schema 为 obs 标准 | 与社区主流对齐；scVI/scANVI/Harmony 全部原生支持 |

## 最近 insight

<!-- 上次 sci-interview / 反思 中提炼的核心洞察，压缩在 5 句以内 -->
<!-- 详细原始对话见 planning/ -->

- 框架本身的科研价值在于"多源 + 多方法并存 + 迭代回跑"，而非单点方法新颖性 → 工程深度比算法新颖更重要
- 与社区主流（scanpy/anndata/CellxGene/scVI）对齐是降低使用门槛与发表成本的关键
- 学生代码主要价值在"具体方法实现细节"（infercnv/Monocle3/UCell 调参），架构层不直接复用
