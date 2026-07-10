---
title: TOM 目录说明
tags: [TOM, cleanup-pending]
created: 2026-07-10
updated: 2026-07-10
---

# TOM — 来源待确认目录

本目录当前仅含一个文件 `scrna_TOM.rda`。

## 现状

截至 2026-07-10,无法从代码或提交历史确认此文件的确切来源与用途。已知线索：

- `TOM` 在 hdWGCNA 分析中指拓扑重叠矩阵（Topological Overlap Matrix）,基因
  共表达模块 notebook `notebooks/07_downstream/15_gene_modules.ipynb` 使用
  hdWGCNA,代码中出现过 `TOM` 变量。
- 但全仓没有任何脚本或 notebook 读取或写入本目录下的 `scrna_TOM.rda` 文件。

## 待办

请项目负责人确认此文件是否为 `15_gene_modules.ipynb` 的某次运行中间产物。若是
一次性中间产物且可重新生成,建议移入 `results/` 并纳入 gitignore,或直接删除
本目录；若有长期保留价值,请在此补充其生成方式与用途。
