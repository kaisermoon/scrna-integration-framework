---
title: 结果目录说明
tags: [results, guide]
created: 2026-07-10
updated: 2026-07-10
---

# results — 分析产出目录

本目录存放流水线运行产生的所有中间与最终结果。**目录下的产出文件全部不进入
git**（见项目根 `.gitignore`）,由各 notebook 运行时生成。

## 子目录约定

| 子目录 | 内容 |
|--------|------|
| `data/` | 各 stage 输出的 h5ad 对象,如 `02_merged.h5ad`、`06_annotated.h5ad`,构成流水线的数据传递链 |
| `figures/` | 图表输出。scanpy 全局图目录已统一指向此处（`sc.settings.figdir`）|
| `tables/` | 表格输出,如差异表达结果、丰度统计等 CSV |

## h5ad 命名与 stage 对应

`data/` 下的 h5ad 文件名以 stage 编号开头,与 `notebooks/` 编号一一对应。上游
stage 的输出是下游 stage 的输入,各 notebook 的具体输入输出契约见
`notebooks/README.md`。
