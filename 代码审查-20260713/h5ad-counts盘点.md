---
title: "01-03 h5ad counts 盘点"
id: h5ad-counts-audit-20260713
tags:
  - project/scRNA-seq整合分析框架
  - code-review
  - data-contract
  - AnnData
created: 2026-07-13
updated: 2026-07-13
---

# 01-03 h5ad counts 盘点

## 结论

> [!important] 核心发现
> Nowicki 当前 manifest 实际读取的 `data/_subset/nowicki/nowicki_subset.h5ad` **并不缺 raw counts**：`adata.raw.X` 是与 `adata.X` 同形、同基因顺序的非负整数矩阵；`adata.X` 则由 `uns["X_name"] = "logcounts"` 明确标记为归一化表达。当前问题是 01 notebook 没有把 `raw.X` 提升为统一的 `layers["counts"]`，而非源文件没有 counts。

当前四个 01 来源均有 counts 入口：

| 来源 | counts 位置 | 当前证据 | 正式准入判断 |
|---|---|---|---|
| Nowicki | 输入 h5ad 的 `raw.X` | **高**：全量值检查 + `X_name=logcounts` + raw/current 基因完全对齐 | 当前 2,500-cell manifest 夹具可准入；完整原文件恢复后仍应按同一规则复核 |
| Kim | 输入 h5ad 的 `X` | **高**：全量值检查为非负整数 | 可准入，01 应复制到 `layers["counts"]` 后再统一生成 normalized `X` |
| Nancang | filtered 10x MTX；另有 raw 10x MTX | **高**：现存 01/02 h5ad 的 `X` 全量为非负整数，且六个样本的 filtered/raw MTX 均存在 | filtered counts 可准入；raw droplet pool 是否满足 SoupX 是另一项独立审查问题 |
| Yue | 29 个 `*_count.txt.gz` | **高**：全量 15,061,300 个表达值均为非负整数 | 可准入，读入后应建立 `layers["counts"]` |

因此，raw-counts 准入制度不要求立即剔除 Nowicki。更准确的整改目标是：**四来源在 01 输出时统一建立真实 `layers["counts"]`，并用统一流程生成 `X`；禁止在 03 从 normalized `X` 反向伪造 counts。**

## 范围与方法

本盘点承接[[全仓代码审查总报告]]，只覆盖 01-03 的实际输入链及现存相关 h5ad，不审查 04 之后的分析产物。

- 定位依据：四个 01 per-dataset notebook、四份 manifest、`02_merged.ipynb`、`03_normalized.ipynb`。
- h5ad 检查：shape、dtype、稀疏性、min/max、非有限值、负值、全量整数性、`raw.X`、layers、var names 对齐以及 `uns` 预处理元数据。
- “整数性”定义：每个显式存储值与最近整数偏差不超过 `1e-6`；下表均为全量检查，不是抽样。
- float32 仅表示存储 dtype；若所有值均精确为非负整数，仍判为 count-compatible。正式契约可另要求整数 dtype，但不能据 float32 单独否定 counts。

## 实际输入链与文件状态

01 当前 manifest 输入：

- Nowicki：`data/_subset/nowicki/nowicki_subset.h5ad`
- Kim：`data/_subset/kim/kim_subset.h5ad`
- Nancang：`data/_subset/nancang/` 下六个样本的 10x MTX，不直接输入 h5ad
- Yue：`data/_subset/yue/` 下 29 个 count 文本，不直接输入 h5ad

02 notebook 预期读取：

- `results/01_nancang_v1.h5ad`
- `results/01_kim_v1.h5ad`
- `results/01_nowicki_v1.h5ad`
- `results/01_yue_v1.h5ad`

03 notebook 预期读取 `results/02_merged_v1.h5ad`。截至 2026-07-13，上述四个正式 01 输出、02 输出及正式 03 输出在当前工作树中均不存在，无法核验其运行后容器状态。

另有一条旧命名 Nancang 01-03 链，以及一个供 03+ 测试的旧夹具：

- `results/nancang_01_loaded_v1.h5ad`
- `results/nancang_02_qcd_v1.h5ad`
- `results/nancang_03_normalized_v1.h5ad`
- `data/_subset/fixture_B_qcd.h5ad`

## h5ad 全量检查结果

| 文件 / 矩阵 | shape | dtype | min-max | 负值 | 非整数值 | 判定 |
|---|---:|---|---:|---:|---:|---|
| Nowicki `X` | 2,500 x 25,853 | sparse float32 | 0-9.5626 | 0 | 4,234,398 / 4,234,398 nnz | normalized/logcounts，不是 counts |
| Nowicki `raw.X` | 2,500 x 25,853 | sparse float32 | 0-37,947 | 0 | 0 / 4,234,398 nnz | raw counts，count-compatible |
| Kim `X` | 1,498 x 38,606 | sparse float32 | 0-2,640 | 0 | 0 / 5,815,711 nnz | raw counts，count-compatible |
| fixture B `X` | 5,998 x 25,593 | sparse float32 | 0-9.0710 | 0 | 10,773,917 / 10,773,921 nnz | normalized/logcounts，不是 counts |
| fixture B `raw.X` | 5,998 x 25,853 | sparse float32 | 0-25,745 | 0 | 0 / 10,773,928 nnz | raw counts，count-compatible |
| Nancang 01 `X` | 1,500 x 38,606 | sparse float32 | 0-24,690 | 0 | 0 / 3,972,583 nnz | raw counts，count-compatible |
| Nancang 02 `X` | 1,091 x 38,606 | sparse float32 | 0-24,690 | 0 | 0 / 2,730,175 nnz | filtered raw counts，count-compatible |
| Nancang 03 `X` | 1,091 x 38,606 | sparse float32 | 0-8.9217 | 0 | 2,730,175 / 2,730,175 nnz | normalized/log1p，不是 counts |
| Nancang 03 `layers["counts"]` | 1,091 x 38,606 | sparse float32 | 0-24,690 | 0 | 0 / 2,730,175 nnz | raw counts，count-compatible |

以上矩阵均未发现 NaN 或 Inf。

## Nowicki 专项证据

### 1. `raw.X` 与 `X` 的语义明确分离

- `uns["X_name"] = "logcounts"`，直接说明当前 `X` 是 log-normalized expression。
- `X` 的所有 4,234,398 个非零值均为非整数。
- `raw.X` 的所有 4,234,398 个非零值均为精确非负整数。
- `X` 与 `raw.X` 的稀疏非零位置完全一致，符合“同一 counts 矩阵经过归一化/对数变换”的结构特征。
- `raw.X` 每细胞总 counts 范围为 284-198,097，中位数 3,788.5，符合 UMI count 的数量级与离散性。

### 2. 基因对齐完整

- 输入时 `raw.var_names` 与 `adata.var_names` 完全相等，均为 25,853 个 Ensembl 基因。
- `raw.X` 和 `X` shape 完全一致，无需做基因交集或补零。
- 01 notebook 随后调用 `sync_gene_ids()` 将当前 `adata.var.index` 从 Ensembl 改为 symbol；该操作不重排矩阵。为避免 raw/current var names 在命名转换后表面不一致，最稳妥的执行点是在读入并验证后立即把 `raw.X.copy()` 放入 `layers["counts"]`，再进行基因名同步。

### 3. 来源元数据

文件带有 CELLxGENE schema 7.0.0 元数据、论文及数据版本 citation，并将 `X` 命名为 `logcounts`。这与数值检查相互印证。证据足以支持当前夹具中的 `raw.X` 作为真实 counts 使用。

### 4. 尚未覆盖的边界

生成夹具时记录的完整 Nowicki 原文件位于外部只读数据目录，但该路径在当前机器上不存在，因此不能直接复核完整 293,823-cell 文件。当前结论严格适用于 manifest 实际读取的 2,500-cell 文件；恢复完整文件后，应复跑同一组 shape、整数性、非负性、raw/current var 对齐检查。

## 其他来源与旧链

### Kim

Kim h5ad 没有 `raw`、layers 或 `uns` 元数据，但 `X` 全量为非负整数 counts。缺的是显式语义标记，不是 counts 本身。01 读入后应先验证并复制为 `layers["counts"]`，再统一归一化。

### Nancang

六个样本均同时存在 `filtered_feature_bc_matrix/` 与 `raw_feature_bc_matrix/`。旧版 01/02 h5ad 的 `X` 保持整数 counts；旧版 03 已正确保留 `layers["counts"]`，同时把 `X` 转为 log1p 表达。

但 raw MTX 是否保留足够的空液滴用于 SoupX，不由“数值是整数”这一检查回答；该问题继续沿用[[全仓代码审查总报告]]中的独立 BLOCK，不应与正式整合所需的 filtered counts 混为一项。

### Yue

29 个压缩 count 文本共检查 15,061,300 个表达值：min=0、max=4,956、负值=0、非整数=0、非数值=0。输入不是 h5ad，但具备生成真实 `layers["counts"]` 的充分条件。

### fixture B

fixture B 也保留 `raw.X` counts。其当前 `X` 只含 25,593 个基因，而 raw 含 25,853 个基因；当前 var names 是 raw var names 的完整子集，交集为 25,593。若继续使用该夹具，应按当前 var names 对 raw 取子集后建立 counts layer，不能直接假设两者同 shape。

## 对正式数据契约的建议

01 每来源输出前应满足：

1. `layers["counts"]` 必须存在，并通过非负、有限、整数性、shape 与 var 顺序检查。
2. `X` 必须按项目统一规则从 `layers["counts"]` 生成；来源自带 normalized `X` 不直接参与 02 合并。
3. `uns["counts_provenance"]` 记录 counts 来自 `X`、`raw.X`、10x filtered MTX 或 count 文本，以及验证结果。
4. `uns["preprocessing_by_source"]` 按来源保存，02 不再把 `preprocessing_done` 做集合并集后代表所有来源。
5. 02 concat 前逐来源硬检查 counts layer；03 只消费真实 counts layer，绝不从当前 `X` 复制或取整构造 counts。
6. counts-dependent 方法仅在上述契约通过时启用；失败来源进入探索轨或停止，而不是静默降级。

> [!warning] 当前仍不能解除 01-03 BLOCK
> counts 可得性问题已经找到可行解，但正式四来源 01 输出和 02/03 输出尚未生成，代码也尚未实施统一契约。解除 BLOCK 需要修复后重新运行，并对四个 01 输出、02 merged、03 normalized 逐文件复验。
