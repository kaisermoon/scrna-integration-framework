---
title: "Per-dataset Notebook 编写规范"
tags:
  - scRNA-seq
  - notebook
  - convention
  - data-ingestion
  - qc
created: "2026-07-25"
updated: "2026-07-25"
---

# Per-dataset Notebook 编写规范

本文档汇总了在多个数据集的实际接入与调试过程中积累的通用技术教训，按主题分节。每条规则包含具体做法与失败现象（为什么），面向非计算机专业的研究人员与学生。

---

## 一、数据读取与格式差异

### 1.1 10X 目录结构：双层优先，扁平后备

**规则**：子目录发现逻辑必须按优先级检查两种结构。

```python
for sub in sorted(data_dir.iterdir()):
    if sub.is_dir() and not sub.name.startswith("."):
        # 优先：CellRanger 标准双层结构
        filtered_dir = sub / "filtered_feature_bc_matrix"
        if filtered_dir.exists():
            sub_dirs.append(filtered_dir)
            continue
        # 后备：扁平结构（matrix.mtx.gz 直接位于子目录下）
        if (sub / "matrix.mtx.gz").exists():
            sub_dirs.append(sub)
```

**sample_id 提取需按层级自适应**：`filtered_feature_bc_matrix` 子目录时，样本标识是上一级目录名；扁平结构时是子目录名本身。

**为什么**：CellRanger 标准输出为 `{sample_id}/filtered_feature_bc_matrix/` 双层结构，但公开数据集的目录组织方式不统一。只检查扁平结构会导致双层结构的样本被识别为 0 个。

### 1.2 `sc.read_10x_h5()` 之后立即 `var_names_make_unique()`

**规则**：每个样本用 `sc.read_10x_h5()` 读取后，紧跟一行去重。

```python
adata_sub = sc.read_10x_h5(h5_file, gex_only=True)
adata_sub.var_names_make_unique()  # 10x h5 基因列表可含重复名
```

**为什么**：10x h5 文件的基因列表可能包含重复基因名（如不同基因组区域的同名基因）。若不去重，后续 `anndata.concat(adatas, join="outer")` 会因 var_names 不唯一而抛出 `InvalidIndexError`。

### 1.3 样本命名须取自 GEO `Sample_title`，不用 FTP 文件名

**规则**：manifest 的 sample 命名必须以 GEO series matrix 中的 `Sample_title` 字段为准。

**为什么**：GEO FTP 服务器的文件名（用于 wget 下载）可能截断完整标题中的关键后缀。曾遇到 FTP 文件名中 "AOIMO" 被截断为 "AO"、"AgOIMO" 被截断为 "AgO" 的情况，导致多个样本被错误归类到错误的疾病分组。这种错误在后续分析中极难发现。

### 1.4 Excel 补充材料：声明 `header_row` 偏移

**规则**：manifest 的 `clinical_metadata` 段支持 `header_row` 字段（默认 0）。

```yaml
clinical_metadata:
  file: "data/Dataset/supplementary.xlsx"
  sheet: "Sheet1"
  header_row: 1  # 跳过大标题行，第二行才是列名
```

**为什么**：部分论文补充材料的 Excel 文件第一行是表格标题（如 "Table S1. Clinical information..."），第二行才是列名。直接 `pd.read_excel()` 会把标题行当列名，导致所有字段映射失败。

### 1.5 Excel 读取后字符串字段统一 `.str.strip()`

**规则**：所有从 Excel 读取的字符串列，映射到 `adata.obs` 后统一去除首尾空格。

```python
adata.obs["disease"] = adata.obs["donor_id"].map(clinical_df["Column"]).str.strip()
```

**为什么**：Excel 单元格中的字符串值可能含尾部空格（如 `"Intestinal "`），不做处理时后续字符串匹配（如 `==` 比较、`groupby`）会因隐性的空格而静默失败，不报错但结果全错。

### 1.6 字段映射后立即检查 NaN 数量

**规则**：每次 `adata.obs[field] = adata.obs[key].map(mapping_dict)` 之后，紧接一行 NaN 检查。

```python
n_nan = adata.obs[field].isna().sum()
if n_nan > 0:
    print(f"[WARNING] {field}: {n_nan}/{len(adata)} 个样本映射失败，请检查 mapping key 是否与原始值精确匹配")
```

**为什么**：mapping key 与原始值的任何微小差异（大小写、空格、特殊字符、拼写错误）都会导致全部 NaN，不报错但下游分析建立在缺失值上。出现意外 NaN 时的第一步排查方向：核对该字段的 `unique()` 值与 mapping 的 key 是否逐字符一致。

---

## 二、基因 ID 与 var_names

### 2.1 raw counts 提取必须在基因 ID 同步之前

**规则**：如果你需要从 `adata.raw.X` 提取原始计数，这一步必须在 `sync_gene_ids`（Ensembl ID → gene symbol 转换）**之前**执行。

```python
# 正确顺序
adata.layers["counts"] = adata.raw.X.copy()  # 此时两侧都是 Ensembl ID
sync_gene_ids(adata)                           # 再将 var_names 转为 gene symbol

# 错误顺序（先转换再提取）会导致 ValueError: 基因名称不匹配
```

**为什么**：`adata.raw.var_names` 是独立存储，不会随 `sync_gene_ids` 自动同步。如果先转换基因名再提取 raw，两侧的 var_names 不一致（一边 gene symbol、一边 Ensembl ID），会触发基因名不匹配错误。

### 2.2 提取 raw counts 后立即 `adata.raw = None`

**规则**：凡是从 `.raw.X` 提取计数到 `layers["counts"]` 后，必须清除 `.raw`。

```python
adata.layers["counts"] = adata.raw.X.copy()
adata.raw = None  # 原始数据已存入 layers，无需保留 .raw
```

**为什么**：scanpy 的 `_check_use_raw()` 检测到 `adata.raw is not None` 时，会自动将 `use_raw` 设为 `True`，导致 `sc.tl.score_genes_cell_cycle` 等函数使用 `adata.raw.var_names`（Ensembl ID）去匹配 gene symbol 基因列表。两侧基因名不一致 → 全部匹配失败 → 评分结果为全 NaN，不报错但数据静默损坏。

### 2.3 细胞周期评分显式 `use_raw=False`

**规则**：调用 `sc.tl.score_genes_cell_cycle` 时必须显式传入 `use_raw=False`。

```python
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, use_raw=False)
```

**为什么**：细胞周期评分基于 Tirosh et al. 2015 的方法，通过归一化的 log 表达量（`adata.X`）计算，而非原始整数计数。即使用了 2.2 规则清除了 `.raw`，显式写出 `use_raw=False` 是防御性编程——万一 `.raw` 被下游代码意外重新赋值，此参数可以避免静默的 NaN。

---

## 三、AnnData 结构与 counts 不变性

### 3.1 `sp.csr_matrix(adata.X, copy=True)` 必须显式指定 `copy=True`

**规则**：所有从 `adata.X` 创建稀疏矩阵副本的代码，必须加上 `copy=True`。

```python
# 正确
adata.layers["counts"] = sp.csr_matrix(adata.X, dtype=np.float32, copy=True)

# 错误：copy 默认 False，layer 与 adata.X 共享内存
adata.layers["counts"] = sp.csr_matrix(adata.X, dtype=np.float32)
```

**为什么**：`sp.csr_matrix` 的默认行为是 `copy=False`，导致 `layers["counts"]` 与 `adata.X` 共享底层的 `indptr` 数组。下游 `sc.pp.calculate_qc_metrics(layer="counts")` 内部调用 `eliminate_zeros()` 会原地修改共享的 `indptr`，同时损坏 `adata.X`，使得 counts checksum 门禁检测到未预期的数据变化。

### 3.2 从 dense array 创建 AnnData 必须先转稀疏

**规则**：从 pandas DataFrame 或 numpy array 创建 AnnData 时，必须先 `sp.csr_matrix()` 包装。

```python
# 正确
adata = anndata.AnnData(X=sp.csr_matrix(df.T.values.astype(np.float32)))

# 错误：X 是 dense numpy array，checkpoint 的稀疏门禁会先于 expression_contract 检查到并报 FAILED
adata = anndata.AnnData(X=df.T.values.astype(np.float32))
```

**为什么**：checkpoint 的 `x_sparse_float32` 门禁在 expression_contract cell 中检查 `sp.issparse(adata.X)`，如果你等到 expression_contract 才转换就太晚了——门禁会先于转换执行。必须在 AnnData 创建时就已是稀疏矩阵。

### 3.3 多样本 `join="outer"` 拼接后必须 `fillna(0)`

**规则**：凡是用 `join="outer"` 拼接多个样本的矩阵，后续必须检查并填充 NaN。

```python
all_cells_df = pd.concat([df1, df2, ...], axis=1, join="outer")
nan_count = all_cells_df.isna().sum().sum()
if nan_count > 0:
    all_cells_df = all_cells_df.fillna(0)
```

**为什么**：不同样本的基因集不尽相同时（如某样本仅 11k 基因，其他 15k），`join="outer"` 用 NaN 填充缺失基因。NaN 进入 X 矩阵后，`calculate_qc_metrics` 对含 NaN 的细胞计算出 NaN 的 `total_counts`，连锁导致 `log_complexity` 全 NaN，最终在 `np.percentile` 中抛 `IndexError`。填充 0 意为"该基因在此样本中未检测到"。

### 3.4 `anndata.concat` 的 `index_unique` 陷阱

**规则**：只有当不同样本的 barcode 本身可能重复（无前缀区分）时，才使用 `index_unique` 参数。

```python
# 正确：obs_names 已由样本前缀保证唯一，不传 index_unique
adata = anndata.concat(adatas, join="outer")

# 错误：传了 index_unique="_"，会追加 _0/_1 批次后缀
adata = anndata.concat(adatas, join="outer", index_unique="_")
```

**为什么**：anndata 0.12+ 版本中，即使 obs_names 已唯一（每个样本有唯一前缀），`index_unique="_"` 仍会对所有 obs_names 追加 `_0`/`_1` 批次后缀。后果是 SoupX 等需要按原始 barcode 匹配的步骤因格式不一致而全部失败（0% 映射率）。

---

## 四、QC 流程与阈值

### 4.1 MT% 自适应阈值必须设置硬上限

**规则**：在 per-sample MAD 自适应阈值计算后，对 MT% 上限做 clamp。

```python
MAX_PCT_MT_HARD = 40  # 推荐值，根据组织类型调整

for sample_id in thresholds["pct_counts_mt"]["per_sample"]:
    original = thresholds["pct_counts_mt"]["per_sample"][sample_id]["upper"]
    clamped = min(original, MAX_PCT_MT_HARD)
    if clamped < original:
        print(f"  [CLAMP] {sample_id}: MT% {original:.1f}% → {clamped:.1f}%")
    thresholds["pct_counts_mt"]["per_sample"][sample_id]["upper"] = clamped
```

**为什么**：MAD 假设数据在 median 周围对称分布。当整个样本的 MT% 都整体偏高（而非长尾），MAD 的中心估计本身已不可靠。N_MAD=5 可推出无意义的 60%+ 上限，几乎不过滤任何细胞，失去 QC 意义。硬上限是 MAD 的安全阀而非替代——n_genes 和 total_counts 分布通常对称，仍用 MAD；MT% 是需要特殊处理的指标。

### 4.2 `calculate_qc_metrics(qc_vars=["mt"])` 的前置条件

**规则**：调用 `calculate_qc_metrics` 之前，`adata.var["mt"]` 列必须已存在。

```python
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], layer="counts", inplace=True)
```

**同时注意列名差异**：`calculate_qc_metrics` 产出的列名是 `n_genes_by_counts`，如果后续代码依赖 `n_genes`，需添加别名。

```python
adata.obs["n_genes"] = adata.obs["n_genes_by_counts"]
```

**为什么**：`qc_vars=["mt"]` 要求 `var` 中已存在对应布尔列，否则 `KeyError`。列名 `n_genes_by_counts`（而非 `n_genes`）是 scanpy API 的标准命名；对于从 h5ad 加载的预计算数据这些列通常已存在，但对于 10X 原始数据需手动创建。

### 4.3 SoupX 应在 doublet 检测之前执行

**规则**：10X 原始数据的处理流程中，SoupX 环境 RNA 校正应紧接在初始 QC 之后、doublet 检测之前。

**为什么**：环境 RNA 会人为抬高某些基因（含线粒体基因）的表达水平，如果不先校正就直接做 doublet 检测，污染会干扰 QC 阈值判断与 doublet 评分的准确性。正确的处理顺序：加载数据 → 初始 QC → SoupX → QC 重算 → 自适应阈值 → doublet 检测 → 细胞周期评分 → QC 过滤。

---

## 五、SoupX 工程

### 5.1 SoupX 直接读原始 10X 目录

**规则**：SoupX R 脚本的 filtered 和 raw 输入必须直接指向 CellRanger 原始输出目录。

```python
# 正确：R 脚本直接读原始目录
filtered_mtx_dir = data_dir / sample_id / "filtered_feature_bc_matrix"
raw_mtx_dir = data_dir / sample_id / "raw_feature_bc_matrix"

# 错误：Python 先导出临时文件再传路径，barcode 格式不一致
```

**为什么**：若 Python 先将 AnnData 写出临时文件再传给 R 脚本，barcode 格式可能因 obs_names 重命名而改变（R 读到 `AAACCCAAGATACCAA-1`，Python 中是 `{sample_id}_AAACCCAAGATACCAA-1`），导致校正后的 barcode 与 adata.obs_names 100% 不匹配。直接从同一 CellRanger 输出目录读取，两侧 barcode 格式天然一致。

---

## 六、checkpoint 与 run 契约

### 6.1 run_id 一律使用 `prepare_run()` 返回的实际值

**规则**：所有引用 run_id 的位置（manifest_payload、adata.uns）必须使用 `prepare_run()` 的返回值，不使用全局变量。

```python
run_paths = prepare_run(...)

# 正确
manifest_payload = {"run_id": run_paths.run_id, ...}
adata.uns["run_id"] = run_paths.run_id

# 错误：使用全局变量 RUN_ID
manifest_payload = {"run_id": RUN_ID, ...}
```

**为什么**：`prepare_run()` 实现了自动递增——若 `run001` 已存在，自动尝试 `run002`。全局变量 `RUN_ID` 仍是硬编码的 `"run001"`，与 `prepare_run()` 实际使用的 run_id 不一致，导致 `promote_run()` 的校验逻辑报错。属于"不报错时数据静默错误、报错时信息误导"的双重陷阱。

### 6.2 `numpy.bool_` 必须先转为 Python `bool` 才能 JSON 序列化

**规则**：checkpoint 的 `hard_postconditions` 构建完成后，统一转换类型。

```python
hard_postconditions = {
    k: bool(v) if isinstance(v, (bool, np.bool_)) else v
    for k, v in hard_postconditions.items()
}
```

**为什么**：`numpy.bool_` 继承自 `numpy.generic`，不是 Python 原生 `bool` 的子类。标准库 `json.dump` 不认，抛 `TypeError: Object of type bool_ is not JSON serializable`。

### 6.3 counts checksum 的语义：QC 过滤区间内的不变性

**规则**：counts checksum 在 expression_contract cell 建立（基于全量 adata）；SoupX 完成后立即断言 counts 未被中间步骤改写（红线门禁）；QC 过滤后更新 checksum 为过滤后的值。

```python
# 过滤完成后更新 checksum
_counts_checksum = float(adata.layers["counts"].sum())
_counts_checksum_nnz = int(adata.layers["counts"].nnz)
```

**为什么**：checksum 比对的是"QC 过滤区间内 counts 未被 SoupX 等中间步骤改写"，而非"全流程 counts 不变"。QC 过滤时 `adata = adata[keep].copy()` 是合法的细胞子集化（counts 自然减少），不是 counts 被篡改。如果不更新 checksum，checkpoint 会报假阳性失败。

---

## 七、manifest 约定

### 7.1 `obs_mapping` 只写需要重命名的字段

**规则**：manifest 的 `obs_mapping` 段只列出原始字段名与目标字段名不同的映射。

```yaml
# 正确
obs_mapping:
  "patient": "donor_id"       # 需要重命名
  "cell_type_original": "cell_type"  # 需要重命名

# 错误：自我映射会触发"复制后删除源列"逻辑，导致字段被删
obs_mapping:
  "disease": "disease"        # 源==目标，会被误删！
```

**为什么**：notebook 中执行 obs_mapping 的代码逻辑是"复制目标列后删除源列"。当 source 和 target 相同时，复制完立即删掉——字段消失，下游验证报告字段"缺失"，根因在于 manifest 写了不该写的自我映射。

---

## 八、其他常见陷阱

### 8.1 `sc.read_10x_h5()` 重复基因名

见 1.2。此条独立列出因其出错率极高：每次 `read_10x_h5` 调用后都需要这行，漏一次就足以触发 `InvalidIndexError`。

### 8.2 从文件名提取元数据的风险

不要从 FTP 文件名或本地文件名提取样本名称、疾病分组等关键信息。文件名的截断、编码、拼写均不在你的控制范围内。所有样本元数据应来自 GEO series matrix 的 `Sample_title` 字段，或论文补充材料的正式表格。

---

## 编写检查清单

新数据集 notebook 编写完成后，逐条确认：

- [ ] 子目录发现同时支持 `filtered_feature_bc_matrix/` 双层和扁平结构
- [ ] `sc.read_10x_h5()` 之后有 `var_names_make_unique()`
- [ ] sample 命名来源为 GEO `Sample_title`（非 FTP 文件名）
- [ ] Excel 读取检查了 `header_row` 偏移
- [ ] 字符串字段已 `.str.strip()`
- [ ] 字段映射后有 NaN 数量检查
- [ ] raw counts 提取在 sync_gene_ids 之前
- [ ] 提取后 `adata.raw = None`
- [ ] 细胞周期评分显式 `use_raw=False`
- [ ] `sp.csr_matrix(..., copy=True)`
- [ ] AnnData 创建时 X 已是稀疏矩阵
- [ ] `join="outer"` 拼接后 `fillna(0)`
- [ ] obs_names 已唯一下不传 `index_unique`
- [ ] MT% 自适应阈值有硬上限 clamp
- [ ] `calculate_qc_metrics(qc_vars=["mt"])` 前有 `var["mt"]` 赋值
- [ ] 注意 `n_genes_by_counts` vs `n_genes` 列名差异
- [ ] SoupX 在 doublet 检测之前
- [ ] SoupX 直接读原始 10X 目录
- [ ] `run_id` 使用 `prepare_run()` 返回值
- [ ] `hard_postconditions` 有 `numpy.bool_` 类型转换
- [ ] QC 过滤后更新 counts checksum
- [ ] `obs_mapping` 无自我映射条目
