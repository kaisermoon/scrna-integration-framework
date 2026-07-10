---
title: 下游分析模块导览
tags: [notebooks, downstream, guide]
created: 2026-07-10
updated: 2026-07-10
---

# 07_downstream — 下游分析模块

本目录下每个 notebook 是一个独立的下游分析模块,全部以 `06_annotated.ipynb`
的输出 `results/data/06_annotated.h5ad` 为输入,彼此之间没有执行先后依赖,可
按需单独运行。上游流水线（01 到 06）的总览见上一级目录的 `notebooks/README.md`。

## 模块清单

| Notebook | 分析内容 | 主要工具 | 运行环境 |
|----------|---------|---------|---------|
| `07_deg.ipynb` | 细胞类型间差异表达基因 | scanpy rank_genes_groups | Python |
| `08_pseudobulk_deg.ipynb` | 样本级拟批量差异表达 | decoupler + PyDESeq2 | Python |
| `09_cnv.ipynb` | 拷贝数变异推断 | infercnvpy | Python |
| `10_pseudotime.ipynb` | 拟时序（扩散拟时 DPT） | scanpy | Python |
| `10b_pseudotime_monocle3.ipynb` | 拟时序（Monocle3） | Monocle3 | R |
| `10c_pseudotime_cellrank2.ipynb` | 拟时序（CellRank2） | cellrank | Python |
| `10d_pseudotime_cytotrace2.ipynb` | 分化潜能打分 | CytoTRACE2 | R |
| `10e_pseudotime_compare.ipynb` | 上述拟时序方法横向比较 | 交叉比对 | Python |
| `11_abundance.ipynb` | 细胞类型丰度组间差异 | scCODA | Python（scrna-sccoda 环境）|
| `12_pathway.ipynb` | 通路富集 | decoupler + MSigDB | Python |
| `13_grn.ipynb` | 基因调控网络 | pySCENIC | Python |
| `14_cell_communication.ipynb` | 细胞间通讯 | LIANA | Python |
| `15_gene_modules.ipynb` | 基因共表达模块 | hdWGCNA | R |
| `16_trajectory_de.ipynb` | 轨迹相关差异表达 | tradeSeq | R |

## 拟时序列的产生与消费（重要,防止踩坑）

拟时序有多个 notebook 各自产生一列 obs,下游 notebook 又要选用其中一列。历史上
多次出现「选用了一个根本没有任何 notebook 产生的列名」这类错误。下表固定当前
各 pseudotime 列的产生方与消费方,改动前务必先核对。

| obs 列名 | 由哪个 notebook 写入 | 状态 |
|----------|---------------------|------|
| `pseudotime_monocle3_v1` | `10_pseudotime.ipynb` / `10b_pseudotime_monocle3.ipynb` | 真实产生 |
| `cellrank2_pseudotime` | `10c_pseudotime_cellrank2.ipynb` | 真实产生 |
| `dpt_pseudotime` | 暂无 notebook 产生 | 保留为占位列名,将来实现 scanpy DPT 时零改动即可启用 |

- `10_pseudotime.ipynb` 的候选列表是「本 notebook 产物的可视化选择器」。
- `10c_pseudotime_cellrank2.ipynb` 的候选列表是「CellRank2 kernel 的输入选择器」,
  在写入 `cellrank2_pseudotime` 之前执行,因此**不包含**它自己稍后才产生的这一列。
- `16_trajectory_de.ipynb` 的候选列表是「下游消费选择器」,包含上面所有真实产生的列。

新增拟时序方法时,先在本表登记新列名与产生方,再在消费方 notebook 的候选列表里
加入该列,避免出现无产生方的空列名。
