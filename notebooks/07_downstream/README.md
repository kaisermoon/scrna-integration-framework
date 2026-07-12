---
title: 下游分析模块导览
tags: [notebooks, downstream, guide]
created: 2026-07-10
updated: 2026-07-12
---

# 07_downstream — 下游分析模块

本目录下每个 notebook 是一个独立的下游分析模块,全部以 `06_annotated.ipynb`
的输出 `results/data/06_annotated.h5ad` 为输入,彼此之间没有执行先后依赖,可
按需单独运行。上游流水线（01 到 06）的总览见上一级目录的 `notebooks/README.md`。

## 命名约定

本目录内文件用 `D` 前缀编号（`D01`、`D02`……）。`D` 表示 downstream,编号
只是清单序位,**不表示执行顺序**——各模块彼此独立、均以 `06_annotated.h5ad`
为唯一输入。用 `D` 前缀而非续接主线的 `07`、`08`……,是为了避免"目录序位
（07_downstream）与目录内文件序位（07_deg）语义重叠"造成的困惑。

## 模块清单

| Notebook | 分析内容 | 主要工具 | 运行环境 |
|----------|---------|---------|---------|
| `D01_deg.ipynb` | 细胞类型间差异表达基因 | scanpy rank_genes_groups | Python |
| `D02_pseudobulk_deg.ipynb` | 样本级拟批量差异表达 | decoupler + PyDESeq2 | Python |
| `D03_cnv.ipynb` | 拷贝数变异推断 | infercnvpy | Python |
| `D04_pseudotime.ipynb` | 拟时序综合（转录组熵 + CytoTRACE + root 识别 + Monocle3） | cellrank + Monocle3(R) | Python + R |
| `D05_pseudotime_monocle3.ipynb` | 拟时序（Monocle3 独立版） | Monocle3 | R |
| `D06_pseudotime_cellrank2.ipynb` | 拟时序（CellRank2 完整 fate-mapping） | cellrank | Python |
| `D07_potency_cytotrace2.ipynb` | 分化潜能打分（非拟时序） | CytoTRACE2 | R |
| `D08_pseudotime_compare.ipynb` | 上述拟时序/潜能方法横向比较 | 交叉比对 | Python |
| `D09_abundance.ipynb` | 细胞类型丰度组间差异 | scCODA | Python（scrna-sccoda 环境）|
| `D10_pathway.ipynb` | 通路富集 | decoupler + MSigDB | Python |
| `D11_grn.ipynb` | 基因调控网络 | pySCENIC | Python |
| `D12_cell_communication.ipynb` | 细胞间通讯 | LIANA | Python |
| `D13_gene_modules.ipynb` | 基因共表达模块 | hdWGCNA | R |
| `D14_trajectory_de.ipynb` | 轨迹相关差异表达 | tradeSeq | R |

`D04` 到 `D08` 是围绕"细胞分化/轨迹"的一组相关分析:`D04`/`D05`/`D06` 各用
不同方法（综合版 / Monocle3 / CellRank2）估计拟时序,`D07` 用 CytoTRACE2 打
分化潜能（potency,非拟时序）,`D08` 把这几种结果横向比较。编号相邻只表示主题
相关,不表示同一方法的变体,也不表示执行先后。

## 拟时序列的产生与消费（重要,防止踩坑）

拟时序有多个 notebook 各自产生一列 obs,下游 notebook 又要选用其中一列。历史上
多次出现「选用了一个根本没有任何 notebook 产生的列名」这类错误。下表固定当前
各 pseudotime 列的产生方与消费方,改动前务必先核对。

| obs 列名 | 由哪个 notebook 写入 | 状态 |
|----------|---------------------|------|
| `pseudotime_monocle3_v1` | `D04_pseudotime.ipynb` / `D05_pseudotime_monocle3.ipynb` | 真实产生（Monocle3 R 桥） |
| `cellrank2_pseudotime` | `D06_pseudotime_cellrank2.ipynb` | 真实产生 |
| `dpt_pseudotime` | 暂无 notebook 产生 | 保留为占位列名,将来实现 scanpy DPT 时零改动即可启用 |

- `D04_pseudotime.ipynb` 用转录组熵、CytoTRACE 与 Monocle3 R 桥做综合拟时序,
  真实写入的拟时序列是 `pseudotime_monocle3_v1`（来自 Monocle3）。它**不调用**
  scanpy 的 `sc.tl.dpt`,因此**不产生** `dpt_pseudotime`;`dpt_pseudotime`
  仅是留给将来 scanpy DPT 实现的占位列名。
- `D04_pseudotime.ipynb` 的候选列表是「本 notebook 产物的可视化选择器」。
- `D06_pseudotime_cellrank2.ipynb` 的候选列表是「CellRank2 kernel 的输入选择器」,
  在写入 `cellrank2_pseudotime` 之前执行,因此**不包含**它自己稍后才产生的这一列。
- `D14_trajectory_de.ipynb` 的候选列表是「下游消费选择器」,包含上面所有真实产生的列。

新增拟时序方法时,先在本表登记新列名与产生方,再在消费方 notebook 的候选列表里
加入该列,避免出现无产生方的空列名。
