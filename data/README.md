---
title: 数据目录说明
tags: [data, guide]
created: 2026-07-10
updated: 2026-07-10
---

# data — 输入数据目录

本目录存放流水线的输入数据。**目录下的数据文件全部不进入 git**（见项目根
`.gitignore`）,因此全新克隆的仓库中本目录是空的,需要按下述约定手动放入数据。

## 子目录约定

| 子目录 | 内容 | 是否进 git |
|--------|------|-----------|
| `raw/` | 外部原始数据,只读。按 `{第一作者}_{年份}_{GEO编号}/` 组织,例如 `Kim_2025_GSE183904/` | 否 |
| `_subset/` | 测试用抽样子集,由 `scripts/make_test_subset.py` 从 `raw/` 生成,供 pytest 与 CI 使用（约 443 MB）| 否 |

## 放入数据的方式

原始数据由研究者手动放入 `raw/`,或从外部只读数据目录建立软链接。流水线的
第一步（`notebooks/01_per_dataset/`）从 `raw/` 读取对应数据集。测试子集不需要
手动准备,运行 `scripts/make_test_subset.py` 即可从 `raw/` 生成到 `_subset/`。
