# Notebooks 导航

scRNA-seq 整合分析框架的完整分析流水线，按数字编号顺序执行。

## 主管线（按顺序执行）

| 编号 | Notebook | 功能 | 输入 | 输出 |
|------|----------|------|------|------|
| 01 | `01_*_loaded.ipynb` (per-dataset) | 数据加载 + QC | raw data (h5/mtx/h5ad) | `results/01_*_qcd_v1.h5ad` |
| 02 | `02_merged.ipynb` | 多数据集合并 + 合并诊断 | 01 各数据集产物 | `results/02_merged_v1.h5ad` |
| 03 | `03_normalized.ipynb` | 归一化 + HVG 选择 | 02 产物 | `results/03_normalized_v1.h5ad` |
| 04 | `04_embedded.ipynb` | PCA/Harmony/scVI 嵌入 | 03 产物 | `results/04_embedded_v1.h5ad` |
| 05 | `05_clustered.ipynb` | 多分辨率 Leiden 聚类 | 04 产物 | `results/05_clustered_v1.h5ad` |
| 06 | `06_annotated.ipynb` | 多方法注释 + LLM 判决 | 05 产物 | `results/06_annotated_v1.h5ad` |

## 辅助管线（按需执行）

| 编号 | Notebook | 功能 | 触发条件 |
|------|----------|------|----------|
| 06b | `06b_per_cluster.ipynb` | 逐簇深度剖析 | 06 完成后 / 06c 完成后 |
| 06c | `06c_subset.ipynb` | Compartment 子集重分析 | 06 完成后需精细注释亚群 |

## 下游分析（07_downstream/）

下游模块采用 `D` 前缀命名，彼此无执行先后依赖，均以 `06_annotated.h5ad` 为输入、可单独运行。`D` 编号仅为清单序位，不表示执行顺序。目录内导览见 `07_downstream/README.md`。

| 编号 | Notebook | 功能 |
|------|----------|------|
| D01 | `D01_deg.ipynb` | Per-cluster 差异表达基因 |
| D02 | `D02_pseudobulk_deg.ipynb` | Pseudobulk DEG (DESeq2) |
| D03 | `D03_cnv.ipynb` | 拷贝数变异推断 |
| D04 | `D04_pseudotime.ipynb` | 轨迹分析（熵 + CytoTRACE + Monocle3） |
| D05 | `D05_pseudotime_monocle3.ipynb` | 拟时序（Monocle3 独立版） |
| D06 | `D06_pseudotime_cellrank2.ipynb` | 拟时序（CellRank2 命运图谱） |
| D07 | `D07_potency_cytotrace2.ipynb` | 分化潜能打分（CytoTRACE2） |
| D08 | `D08_pseudotime_compare.ipynb` | 拟时序/潜能方法横向比较 |
| D09 | `D09_abundance.ipynb` | 细胞类型丰度比较 |
| D10 | `D10_pathway.ipynb` | 通路富集分析 |
| D11 | `D11_grn.ipynb` | 基因调控网络 |
| D12 | `D12_cell_communication.ipynb` | 细胞通讯 |
| D13 | `D13_gene_modules.ipynb` | 基因共表达模块 |
| D14 | `D14_trajectory_de.ipynb` | Trajectory-关联 DEG (GAM) |

## 典型工作流

```
全局分析:  01(x5) -> 02 -> 03 -> 04 -> 05 -> 06 -> 06b -> D01-D14(按需)

上皮精细分析:  06 -> 06c(epithelial subset) -> 06b(subset mode) -> D04 -> D14

回退重跑:  发现 04 过校正 -> 改 04 PARAMS -> 重跑 05 -> 06 -> ...
          (用 scripts/trace_downstream.py 查看影响链)
```

## 参数调整指南

每个 notebook 的 `PARAMS` cell 集中定义所有可调参数。关键参数跨 stage 传递关系：

- 04 `EMBEDDING_METHODS` -> 决定 04 产出哪些 obsm key
- 04 决策摘要输出 `USE_REP = "..."` -> 复制到 05 `USE_REP`
- 05 `recommended_resolution` (uns) -> 参考后设 06 `LEIDEN_COL`
- 06 `cell_type_final_v1` -> 06b `LABEL_COL` / 06c `SUBSET_FILTER`
