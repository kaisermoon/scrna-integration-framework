---
title: "标记物库使用说明"
type: marker-library-readme
project_id: "scrna-integration-framework"
created: "2026-06-08"
updated: "2026-06-08"
---

# 标记物库（Marker Library）

> **重要声明：本目录下的模板 CSV 不含任何真实标记物基因或文献 PMID。所有标记物内容由 PI（钟子劭，消化科专家）亲自填写。框架只提供格式规范与加载工具。**

## 目录结构

```
references/markers/
├── README.md                              # 本文件
├── _TEMPLATE_markers.csv                  # 通用空模板（仅表头）
├── _TEMPLATE_gastric_epithelial.csv       # 胃上皮标记物模板（仅表头）
└── gastric_epithelial.csv                 # 胃上皮标记物（PI 填写后创建）
    gastric_immune.csv                     # 胃免疫标记物（PI 填写后创建）
    ...                                    # 按 {组织}_{用途}.csv 命名
```

模板文件（`_TEMPLATE_*.csv`）只有表头行，不含任何数据。PI 复制模板文件，改名（去掉 `_TEMPLATE_` 前缀），然后按本 README 的规范逐行填入真实标记物。

## CSV 字段说明

`load_markers()` 读取的 CSV 必须包含以下 6 列（顺序不重要，但列名必须完全一致）：

| 列名 | 类型 | 必填 | 说明 |
|------|------|------|------|
| `tissue` | 文本 | 是 | 组织来源，如 `gastric`、`intestinal`、`liver`。同一文件内可混入多种组织（通常一种组织一个文件）。 |
| `cell_type` | 文本 | 是 | 细胞类型名称。建议使用学界通用缩写，如 `pit_cell`（胃小凹细胞）、`chief_cell`（主细胞）、`SPEM`（解痉多肽表达化生）。同一细胞类型可有多行（每行一个标记物基因）。 |
| `marker` | 文本 | 是 | **基因符号（Gene Symbol）**，如 `MUC5AC`、`TFF2`。必须为 HUGO Gene Nomenclature Committee (HGNC) 官方符号。**注意：不同数据集的基因集不同，使用前必须检查基因是否存在于 `adata.var_names` 中**（见下方「基因存在性检查」节）。 |
| `role` | 枚举 | 是 | 标记物在细胞类型鉴定中的角色。取值必须为以下三者之一：`canonical`（经典标记物，用于阳性鉴定）、`optional`（辅助标记物，可增强置信度但非必需）、`negative`（阴性标记物，该细胞类型应**不表达**此基因，用于排除） |
| `reference` | 文本 | 是 | 文献来源。格式：`<第一作者姓氏> <年份>` 或 **`PMID:XXXXXXXX`（推荐）**。**硬要求：每个 PMID 必须可在 PubMed 核实**（https://pubmed.ncbi.nlm.nih.gov/XXXXXXXX/）。禁止填入无法核实的引用（如"内部数据""个人通讯""实验室经验"）。可接受"综述引用 + 原始文献 PMID"双标注。 |
| `notes` | 文本 | 否 | 备注（可选）。可记录：表达特异性强度、适用条件（如"仅在人胃窦表达"）、与其他标记物的互斥关系、是否仅适用于特定 dataset 等。 |

### `role` 字段取值规则（重要）

- **`canonical`**：该细胞类型公认的标志性基因。在 dotplot、基因集评分中作为阳性证据使用。`load_markers()` 默认返回 `canonical` + `optional`。
- **`optional`**：该细胞类型可表达的基因，但不是唯一标志。可用于辅助判断，不适合单独作为分群依据。默认随 `canonical` 一起返回。
- **`negative`**：该细胞类型**不应表达**的基因。用于反向验证——如果某个 leiden 簇高表达某细胞类型的阴性标记物，则排除该细胞类型归属。`load_markers()` 默认**不**返回阴性标记物，需显式指定 `roles=("negative",)`。

**示例（以下为占位符，非真实基因——请删除后填写）**：

```csv
tissue,cell_type,marker,role,reference,notes
gastric,pit_cell,<GENE_SYMBOL>,canonical,PMID:XXXXXXXX,"胃小凹细胞经典标记物"
gastric,pit_cell,<GENE_SYMBOL>,optional,PMID:XXXXXXXX,"辅助标记，非特异性"
gastric,SPEM,<GENE_SYMBOL>,canonical,PMID:XXXXXXXX,"SPEM 核心标记物"
gastric,SPEM,<GENE_SYMBOL>,negative,PMID:XXXXXXXX,"SPEM 不应表达此基因"
```

## 如何填写标记物表（面向消化科专家的操作指南）

### 第 1 步：复制模板

```bash
cp references/markers/_TEMPLATE_markers.csv references/markers/gastric_epithelial.csv
```

### 第 2 步：用 Excel / Numbers / VS Code 打开 CSV

CSV 是纯文本表格格式，任何表格软件均可打开。**推荐 VS Code**（安装 Rainbow CSV 插件后可高亮列），或直接用 Excel。

### 第 3 步：逐行填入标记物

以"胃小凹细胞（pit_cell）"为例，填入流程：

1. 翻查文献（PubMed 检索 `"pit cell" OR "foveolar cell" marker gastric`）
2. 确认该基因确实是 pit_cell 的公认标记物（建议交叉验证 ≥2 篇独立文献）
3. 判断 `role`：
   - 该基因几乎是 pit_cell 的代名词 → `canonical`
   - 该基因在 pit_cell 表达但也在其他上皮细胞表达 → `optional`
   - 该基因在 pit_cell **不**表达，可用于排除 → `negative`
4. 填写 `reference`：**优先填 PMID**（如 `PMID:37995690`），方便后续核查
5. 填写 `notes`（如有需要）：如"仅在胃体小凹细胞高表达"

### 第 4 步：保存，在 notebook 中使用

```python
from scrna_integration import load_markers

# 加载胃上皮标记物（默认返回 canonical + optional）
markers = load_markers("references/markers/gastric_epithelial.csv")

# 在 dotplot 中使用
sc.pl.dotplot(adata, var_names=markers, groupby="leiden")

# 仅在需要阴性标记物时
negative = load_markers("references/markers/gastric_epithelial.csv", roles=("negative",))
```

### 填写注意事项

1. **基因符号准确性**：务必核对 HUGO 官方符号。例如确认是 `MUC5AC` 而非 `MUC5AC`（拼写变体）、是 `TFF2` 而非 `TFF2`（已弃用符号）。
2. **PMID 可核实**：框架不做自动 PMID 验证（避免网络依赖），但 reviewer / PI 会抽查。请确保每个 PMID 对应的论文确实讨论了该基因在该细胞类型中的表达。
3. **一行一个基因**：如果一个细胞类型有 10 个 canonical 标记物，就写 10 行（同一 `cell_type`，不同 `marker`）。`load_markers()` 会自动按细胞类型分组。
4. **不编造**：如果你不确定某个基因是否适合作为标记物，宁可不填。标记物库的质量直接影响 stage 6 注释和下游分析的准确性。

## 基因存在性检查（使用前必须执行）

**不同 scRNA-seq 数据集包含的基因集合不同。** 即使某个基因在文献中是公认标记物，它也可能不在当前数据集的 `adata.var_names` 中。**在使用标记物列表前，必须检查基因是否存在，并报告缺失情况。** 这是 student-code 中反复出现的真实需求（参见 ADR-0008）。

### 标准写法

```python
from scrna_integration import load_markers

markers = load_markers("references/markers/gastric_epithelial.csv")

# 对每个细胞类型，检查其标记物在当前数据集中的存在性
for cell_type, gene_list in markers.items():
    genes_present = [g for g in gene_list if g in adata.var_names]
    genes_missing = [g for g in gene_list if g not in adata.var_names]

    if genes_missing:
        print(
            f"[{cell_type}] {len(genes_missing)}/{len(gene_list)} markers "
            f"not found in dataset: {genes_missing}"
        )

    # 只使用存在的基因做后续分析
    markers[cell_type] = genes_present
```

### 为什么必须检查

- **不同参考基因组**：GRCh38 vs GRCh37 基因符号有差异（如部分基因在 GRCh38 中改名或拆分）。
- **不同测序平台**：10x Genomics 3' vs 5' vs Smart-seq2 捕获的基因集合不同。
- **不同物种**：人类 vs 小鼠基因符号不同（虽然框架已通过 manifest 强制 `species: human`，但标记物库将来可能扩展到小鼠模型）。
- **数据集质量**：低质量细胞或特定组织可能缺失某些基因的表达。

**原则：永远不要假设标记物列表中的所有基因都在当前数据集中存在。先检查，再使用。**

## 加载工具

`load_markers()` 是框架的三大公开函数之一（参见 [ADR-0005](../../docs/adr/0005-load-markers-helper.md) / [SPEC.md](../../SPEC.md)），位于 `src/scrna_integration/markers.py`。

```python
from scrna_integration import load_markers

# 扁平模式（默认）：{cell_type: [gene_list]}
markers = load_markers("references/markers/gastric_epithelial.csv")
# → {"pit_cell": ["MUC5AC", "TFF1", ...], "SPEM": ["TFF2", "MUC6", ...]}

# 仅阴性标记物
negative = load_markers("references/markers/gastric_epithelial.csv", roles=("negative",))

# 完整三层字典：{cell_type: {canonical: [...], optional: [...], negative: [...]}}
full = load_markers("references/markers/gastric_epithelial.csv", roles=None)
```

## 相关文档

- [ADR-0005: `load_markers` as a legitimate third framework function](../../docs/adr/0005-load-markers-helper.md)
- [ADR-0008: 吸收 student-code 下游技术（含基因存在性检查约定）](../../docs/adr/0008-absorb-student-code-by-rewriting.md)
- [SPEC.md § `load_markers`](../../SPEC.md)
- [框架根 README](../../README.md)
