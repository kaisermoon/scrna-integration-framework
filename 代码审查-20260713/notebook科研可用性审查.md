---
title: "scRNA-seq 整合分析框架：Notebook 科研可用性代码审查"
tags: [代码审查, scRNA-seq, notebook, 科研可复现性, 生物医学用户]
created: 2026-07-13
updated: 2026-07-13
---

# Notebook 科研可用性代码审查

## 审查结论

**Verdict：BLOCK（不建议把当前整套 notebook 直接作为正式科研分析模板交付学生或据此产出论文结果）。**

当前版本已经做对了若干重要事情：大多数 notebook 有集中 `PARAMS` 区、中文注释会解释部分方法选择的 why、核心 stage 有检查点和稀疏矩阵断言、许多可选重型方法有环境守卫，04–06 也尝试把方法比较和 PI 判断留在 notebook 中。这些基础值得保留。

但本轮全量静态审查确认存在会直接改变科研结论的高风险问题，尤其是：不同预处理尺度的数据被直接合并、pseudobulk 把不同细胞类型混进同一个 DESeq2 模型、GRN 把所有 HVG 当作 TF、CellChat 使用原始 counts、trajectory-DE 用 R2 阈值冒充显著性，以及 06 自动把 LLM 高置信标签写成最终标签。上述问题不能靠“100 个测试通过”抵消，因为现有测试没有覆盖这些真实生物学分析契约。

### 严重度统计

| 严重度 | 数量 | 含义 |
|---|---:|---|
| Critical | 9 | 可能产生错误科研结论、泄露临床信息或绕过 PI 判断；正式分析前必须修复 |
| Important | 14 | 会显著妨碍调参、调试、复现或造成错误信心；应在交付生物医学用户前修复 |
| Minor | 6 | 文档、命名和可读性问题；可随统一治理处理 |

## 审查范围与覆盖率

- 已建立并检查完整清单：**29/29 个 `.ipynb`（100%）**。
- 现行与辅助 notebook：27 个；deprecated notebook：2 个。
- 同时检查 notebook 相关文件：`notebooks/README.md`、`notebooks/07_downstream/README.md`、`notebooks/_llm_proposer.py`、`scripts/smoke_run_notebooks.py`、`scripts/deseq2_contrast.R`。
- 对所有 29 个 notebook 做了逐 cell 结构检查和 Python AST 静态语法检查；结果为 **29 个文件、0 个语法错误**。
- 未执行全量真实数据分析。scVI、Monocle3、CellRank2、CytoTRACE2、CellChat、hdWGCNA 等依赖真实数据/重型环境的方法，本报告对其代码路径、输入契约、失败语义和科研方法正确性做静态审查；运行时兼容性仍需独立实测。

### Notebook 清单

| 区域 | 已审文件 |
|---|---|
| 设计期/辅助 | `00_propose_obs_manifest.ipynb`、`01_template_10x.ipynb`、`06b_per_cluster.ipynb`、`06c_subset.ipynb` |
| 来源预处理 | `01_kim.ipynb`、`01_nancang.ipynb`、`01_nowicki.ipynb`、`01_yue.ipynb` |
| 核心主管线 | `02_merged.ipynb`、`03_normalized.ipynb`、`04_embedded.ipynb`、`05_clustered.ipynb`、`06_annotated.ipynb` |
| 下游 D01–D14 | `D01_deg.ipynb`、`D02_pseudobulk_deg.ipynb`、`D03_cnv.ipynb`、`D04_pseudotime.ipynb`、`D05_pseudotime_monocle3.ipynb`、`D06_pseudotime_cellrank2.ipynb`、`D07_potency_cytotrace2.ipynb`、`D08_pseudotime_compare.ipynb`、`D09_abundance.ipynb`、`D10_pathway.ipynb`、`D11_grn.ipynb`、`D12_cell_communication.ipynb`、`D13_gene_modules.ipynb`、`D14_trajectory_de.ipynb` |
| Deprecated | `_deprecated/01_loaded.ipynb`、`_deprecated/02_qcd.ipynb` |

## Critical Issues

### C1. 不同数据集处于 raw counts 与 normalized 两种尺度，却被直接拼接成一个矩阵

- **位置**：`notebooks/01_per_dataset/01_nowicki.ipynb:cell 3/11/12`；`notebooks/02_merged.ipynb:cell 4/10/26`；`notebooks/03_normalized.ipynb:cell 3/5/7`
- **问题**：Nowicki 明确是“作者已 QC + normalized”，其 `adata.X` 直接进入 02；其他来源的 `X` 是原始 counts。02 用 `anndata.concat` 直接拼接两种尺度。随后又把各数据集 `preprocessing_done` 做集合并集，只要任一来源含 `normalization`，03 就对整个 merged 对象跳过归一化。03 cell 5 还会无条件把当前混合尺度 `X` 复制成 `layers['counts']`。
- **为何重要**：这会让不同数据集的表达值失去可比性，并把 normalized 值伪装成 counts 交给 scVI、DESeq2、CytoTRACE、CellChat 等依赖原始计数的方法。批次结构、HVG、PCA、聚类及所有下游结果均可能被系统性扭曲。
- **具体改法**：禁止在 02 合并不同表达尺度。每个来源必须显式产出两个契约层：`layers['counts']`（真实整数 counts；无则标记不可用）和 `X`（统一 log-normalized 表达）。如果某公开 h5ad 没有 raw counts，应选择：从原始矩阵重建 counts；或将该数据集限制为只参与允许 normalized 输入的方法；不能把 normalized `X` 当 counts。`preprocessing_done` 必须按来源保存为映射，不能做全局并集后控制全体行为。

### C2. D02 pseudobulk 把不同 cell type 的样本混进同一个 DESeq2 模型

- **位置**：`notebooks/07_downstream/D02_pseudobulk_deg.ipynb:cell 5/6/8`；`scripts/deseq2_contrast.R:34-87`
- **问题**：`dc.pp.pseudobulk(... groups_col=CELL_TYPE_COL)` 产生 `(sample, cell_type)` 级 pseudobulk，但导出的全部 cell type 列一起进入 R；R 侧设计只有 `~ disease`，既不逐 cell type 分析，也不把 cell type 纳入模型。
- **为何重要**：不同细胞类型的巨大表达差异会成为主导变异，疾病效应与细胞类型效应混杂。输出表不能解释为任何一个 cell type 的疾病 DEG，是错误的统计模型。
- **具体改法**：以 cell type 为外层显式循环；每个 cell type 单独建立 genes x biological samples 的矩阵，检查每组独立 donor 数，再运行 `design=~ covariates + disease`。每个结果文件名、火山图和 `uns` 键必须包含 cell type。若计划联合建模，至少使用配对/阻断设计并明确 interaction，但不应作为默认入门路径。

### C3. 06 自动接受 LLM 标签，并把“未标注”伪装成最终标注

- **位置**：`notebooks/06_annotated.ipynb:cell 27/29/30/33/34`
- **问题**：LLM 返回 `confidence=HIGH` 时自动预填并直接写入 `cell_type_final_v1`；未填写任何标签时，用 `Cluster_N` 填满同一个 final 列。随后元数据固定写入 `method_basis="PI manual review..."` 和 `rationale="PI reviewed each cluster..."`，即使 PI 实际没有确认。
- **为何重要**：这绕过了项目明确的“PI 最终判断权”。下游 notebook 只检查 final 列是否非空，会把 LLM 自动标签或 `Cluster_N` 当成真实细胞类型继续做 CNV、丰度、通讯和轨迹分析，形成不可见的证据污染。
- **具体改法**：增加显式 `PI_CONFIRMED = False` 闸门；默认只写 `cell_type_suggested_*`，绝不写 final。只有 PI 在 CSV 中逐簇确认、所有簇无空值且 `PI_CONFIRMED=True` 时才能生成 `cell_type_final_*`。未确认时应阻断正式 checkpoint，或写 `status="draft_unconfirmed"` 到单独草稿文件。元数据必须记录确认人/确认时间/决策 CSV hash，而不是固定声称已审阅。

### C4. D11 把所有 HVG 当作候选 TF，得到的不是可信 GRN

- **位置**：`notebooks/07_downstream/D11_grn.ipynb:cell 7`
- **问题**：`tf_names=_hvg_genes`，即所有 HVG 都被声明为转录因子；输入优先使用 raw counts，未控制 library size。后续把这些链接称为 TF-target 调控关系。
- **为何重要**：大量受体、结构蛋白、代谢酶会被错误当成 TF；raw library-size 信号可制造广泛伪相关。cisTarget 不能替代正确 TF universe 和表达预处理。最终 regulon、hub TF 和机制解释均可能失真。
- **具体改法**：在 PARAMS 暴露并加载物种匹配的 curated TF list（如 pySCENIC 官方 human TF list），打印匹配率和缺失率；GRNBoost2 使用合适的 normalized expression；保留非 HVG 的已知 TF/靶基因所需基因空间。若 TF list 缺失应硬阻断正式 GRN，而不是回退为“all HVG”。

### C5. D12 CellChat 使用 raw counts，且给所有标签无条件加 `C` 前缀

- **位置**：`notebooks/07_downstream/D12_cell_communication.ipynb:cell 8/10`
- **问题**：notebook 声称 CellChat 需要原始 counts，并把 counts 直接传给 `createCellChat`；R 脚本未做 `normalizeData`。同时 `paste0("C", cell_type)` 对所有标签加前缀，而非仅处理纯数字标签。
- **为何重要**：CellChat 的表达过量与通讯概率通常基于规范化表达；raw library size 会影响群体间比较。标签会从 `B cell` 变成 `CB cell`，从 `epithelial` 变成 `Cepithelial`，破坏生物医学可读性和与其他结果的对齐。
- **具体改法**：明确采用 CellChat 官方推荐的 normalized data 路径，并在 notebook 中展示归一化步骤和输入分布诊断；仅对 `grepl("^[0-9]+$", label)` 的标签加安全前缀，同时保存双向 label map 并在输出读回时还原原标签。

### C6. D14 不是显著性检验，却以“DEG/FDR/显著”措辞和阈值进入富集

- **位置**：`notebooks/07_downstream/D14_trajectory_de.ipynb:cell 0/1/3/4/6`
- **问题**：实际方法是 `UnivariateSpline` 的训练内 R2 排序；没有 null model、p-value、FDR、样本/批次效应或 donor 层面重复。`R2_THRESHOLD=0.05` 被解释为“显著”，然后据此做通路富集。PARAMS 注释甚至出现不存在的 `FDR_THRESHOLD`。
- **为何重要**：R2 是拟合优度，不是显著性；在同一批细胞上拟合和评估会乐观偏倚，平滑参数、零膨胀、细胞密度和表达量都会影响排名。称其为 trajectory DEG 会误导论文结果。
- **具体改法**：正式方案使用 tradeSeq、GAM/pyGAM 或带样本效应的模型，输出统计量、p-value、FDR 和 effect size；至少加入 permutation/null、交叉验证 R2、有限值过滤、分支特异分析。若暂保留当前实现，必须改名为“exploratory trajectory association score”，删除 DEG/FDR/显著措辞，禁止直接作为富集的正式输入。

### C7. 00 会把 obs/临床表实际值发送给外部 LLM，却没有隐私闸门

- **位置**：`notebooks/00_propose_obs_manifest.ipynb:cell 0/5/7/10`；`notebooks/_llm_proposer.py:74-107`
- **问题**：prompt 包含 obs 前 5 个非空值和临床表头部实际值，随后直接发送到配置的 Anthropic/OpenAI-compatible provider。没有 PHI/PII 检测、脱敏预览、provider 数据边界提示或 PI 二次确认。
- **为何重要**：真实临床表可能含患者 ID、日期、住院号或罕见组合信息。这是数据隐私红线，而“只取前几行”并不等于脱敏。
- **具体改法**：默认只发送列名、dtype、唯一值计数和经脱敏/枚举化后的类别样本；加入 `ALLOW_EXTERNAL_LLM=False`、敏感列正则扫描、待发送 prompt 完整预览和 `PI_CONFIRMED_NO_PHI=True` 双闸门。临床值默认不发送；本体映射优先走本地规则/本地模型。

### C8. 04 在 counts layer 无效、raw 有效时仍选择无效 layer 给 scVI

- **位置**：`notebooks/04_embedded.ipynb:cell 4/18`
- **问题**：当 `layers['counts']` 存在但非整数，而 `raw.X` 是整数时，代码明确“不覆盖”错误 layer，却仍设置 `_counts_key='counts'`；scVI 最终继续读已判定无效的 layer。
- **为何重要**：scVI 负二项模型要求 counts。该 fallback 给用户显示“发现 raw 整数”，实际却训练在坏数据上，诊断输出与真实行为相反。
- **具体改法**：把数据来源与实际模型输入绑定为同一个对象；若 raw 有效，显式对齐并生成新的只读 `counts_for_scvi` layer，或直接阻断让用户修复。全量/分层抽样验证非负、近整数、无 NaN，并在训练前再次断言。

### C9. SoupX 后不重算 QC 指标，过滤使用校正前的陈旧值

- **位置**：`notebooks/01_per_dataset/01_kim.ipynb:cell 20/24`、`01_nancang.ipynb:cell 20/24`、`01_yue.ipynb:cell 20/24`
- **问题**：SoupX 覆盖 `adata.X` 后，`n_genes`、`total_counts`、`pct_counts_mt` 等 QC 列没有重算；后续 MAD 过滤仍使用校正前指标。双细胞也只标记、不提供显式移除/保留决策参数。
- **为何重要**：校正改变表达矩阵后，过滤标准和实际矩阵不再一致；用户看到的前后 QC 与最终数据可能不匹配。双细胞污染会继续进入聚类，且用户难以发现这是设计选择还是遗漏。
- **具体改法**：明确推荐顺序并保持指标一致：原始 QC/过滤决策、SoupX、重新计算校正后 QC 并展示差异；增加 `REMOVE_PREDICTED_DOUBLETS` 参数、每样本阈值表和移除前后 UMAP/QC 诊断。是否移除必须由 PI 显式选择并记录。

## Important Issues

### I1. 分析逻辑仍隐藏在 `src`，违反 notebook 边界铁律

- **位置**：`01_nancang.ipynb:cell 2/3`、`01_nowicki.ipynb:cell 2/3`、`01_yue.ipynb:cell 2/3` 调用 `scrna_integration.io.sync_gene_ids`；`D03_cnv.ipynb:cell 5` 调用 `inject_genomic_positions`。
- **影响**：基因 ID 映射、冲突处理、基因组位置注入都是会改变基因空间和结果的科研处理，不是纯技术管道。生物医学用户无法在 notebook 中看到 mapping 版本、未映射/一对多冲突、丢失基因及修改规则。
- **改法**：把映射表来源、版本、冲突规则、覆盖率统计和实际转换代码展开到 notebook；可保留“RDS 读取一步”在技术桥，但科研映射判断必须可见、可改、可审计。

### I2. `obs_mapping` 方向说明自相矛盾，可能让 LLM 与 PI 写反映射

- **位置**：`00_propose_obs_manifest.ipynb:cell 13`；`notebooks/_llm_proposer.py:129-166/264-307`
- **影响**：文字写“源列 -> 规范字段”，示例与合并实现却按“规范字段 -> 源列”。新数据接入最容易在此产生静默错列。
- **改法**：全仓只保留一个方向并用真实例子展示；写回前校验 key 属于允许规范字段、value 必须存在于 obs 列，并打印 `source -> target` 的人类可读预览。

### I3. D09 的“多组 Kruskal-Wallis + pairwise”没有真正执行统计检验

- **位置**：`notebooks/07_downstream/D09_abundance.ipynb:cell 10`
- **影响**：代码只计算每个疾病组总体细胞比例范围，未使用 per-sample replicates，也未调用已导入的 `kruskal`/`mannwhitneyu`；却用 `range_pct > 5` 标记“显著差异簇”。同时读取未在 PARAMS 定义的 `DISEASE_COL`，而不是 `CONDITION_COL`。
- **改法**：删除伪统计标题和“显著”措辞；基于现成 `proportion_df` 做真正的 sample-level Kruskal-Wallis、事后 pairwise 和多重校正，或只保留为明确的描述性 effect-size 排序。

### I4. 下游方法失败后普遍仍写同名 checkpoint，容易被误认为分析成功

- **位置**：`D02:cell 13`、`D05:cell 12/13`、`D06:cell 14/15`、`D07:cell 10/11`、`D09:cell 18/19`、`D10:cell 20/22`、`D11:cell 15/17`、`D12:cell 21/24`、`D13:cell 19/22`
- **影响**：只要上游 h5ad 可读，即使核心方法没运行，也会生成看似成功的下游文件。文件存在检查无法区分“成功结果”与“原样转存 + skip 标记”。
- **改法**：定义四态：`success / skipped_by_user / unavailable / failed`。只有 success 才写正式 `OUTPUT_PATH`；其他状态写诊断 JSON 或 `_incomplete.h5ad`。结尾用醒目的验收表列出核心结果字段是否存在，失败状态返回硬错误或明确停止下游消费。

### I5. D07 明确写着“API 最佳猜测、未实跑”，不应列为可用正式 stage

- **位置**：`D07_potency_cytotrace2.ipynb:cell 6/7/8/11`
- **影响**：权重检测只要目录非空就判定 ready；API 和返回字段靠多种猜测匹配。用户可能花费大量算力后才发现接口不匹配，或错误映射结果列。
- **改法**：在支持版本环境做一次真实最小数据集验证并锁定 API；加版本断言和官方输出 schema 校验。未验证前将 notebook 标记 `experimental/unvalidated`，从“23-stage 已打通”的正式清单中剥离。

### I6. D06 声称会产生 `cellrank2_pseudotime`，但没有硬输出契约

- **位置**：`D06_pseudotime_cellrank2.ipynb:cell 11`；`notebooks/07_downstream/README.md:51-64`
- **影响**：只有当 `GPCCA` 对象恰好存在 `pseudotime` 属性才写列，异常被空 `except` 吞掉；README 却把该列列为真实产生。D08/D14 可能因此静默跳过 CellRank 方法。
- **改法**：明确 CellRank notebook 的核心产物是 fate probabilities/terminal states，若需要 pseudotime，使用有明确语义的方法生成并验证；结尾断言文档声明的字段实际存在，或同步修改 README。

### I7. D01 跨疾病 DEG 仍在 cell 层面比较且混合细胞类型

- **位置**：`D01_deg.ipynb:cell 9/10`
- **影响**：intro 虽承认 pseudoreplication，但疾病对比仍对所有细胞直接 Wilcoxon，没有按 cell type 分层，也没有 donor 作为独立重复。细胞组成差异会被误解释为基因表达差异。
- **改法**：把该模块明确降级为 exploratory visualization；正式疾病 DEG 路由到修正后的 D02，并按 cell type + donor 做 pseudobulk。若保留 cell-level，至少按 cell type 分层且不报告推断性 p-value。

### I8. D10 ORA 背景和输入基因集不透明，通路活性使用 raw counts

- **位置**：`D10_pathway.ipynb:cell 6/10/16`
- **影响**：ORA 直接取 top N 基因，不按效应量/显著性筛选，也未使用“本数据中可检测基因”作为背景；decoupler 优先 raw counts，容易受 library size 影响。
- **改法**：PARAMS 暴露 DEG 阈值、方向、背景 universe；打印每簇输入基因数和背景数。通路活性使用明确的 log-normalized matrix，并在 sample/cell type 层面聚合后比较，避免把单细胞平均当独立重复。

### I9. D13 Python 层次聚类的距离输入与 Ward 方法不匹配

- **位置**：`D13_gene_modules.ipynb:cell 7`
- **影响**：`sch.linkage(1 - abs(rho_matrix), method='ward')` 把方阵当作 observation matrix，而非 condensed distance；且 Ward 要求欧氏距离，`1-|rho|` 不满足该假设。展示的“朴素模块”树可能不对应预期相关结构。
- **改法**：用 `squareform` 传 condensed distance，并选择适合相关距离的 average/complete linkage；或删掉该临时模块划分，只保留正式 hdWGCNA 结果。

### I10. 04 暴露的部分参数没有真正控制实现

- **位置**：`04_embedded.ipynb:cell 2/18/24`
- **影响**：`SCVI_GENE_LIKELIHOOD` 定义后模型仍硬编码 `zinb`；`SCVI_EARLY_STOPPING` 定义后训练仍硬编码 True；scVI/torch 全局随机性未完整设置。用户修改参数后可能以为结果已变化，实际没有。
- **改法**：每个 PARAMS 变量必须有静态消费检查；运行开始打印“最终生效参数”；设置 `scvi.settings.seed`/相关后端 seed，并记录确定性限制。

### I11. 多处大 cell 不利于断点调试和重跑

- **位置**：`06b:cell 7`（306 行）、`D13:cell 13`（241 行）、`D05:cell 9`（204 行）、`D04:cell 20`（197 行）、`06_annotated:cell 3`（164 行）、`D12:cell 10`（184 行）、`06c:cell 15`（184 行）等。
- **影响**：一个步骤中同时做准备、运行、解析、绘图和写盘；失败只能整 cell 重跑，难以定位中间状态，也不适合医学用户逐步验证。
- **改法**：按“准备输入 -> 显示输入诊断 -> 执行方法 -> 读取结果 -> QC 图 -> PI 决策 -> 保存”拆成 20–60 行 cell。R 脚本可保留可见，但单独放代码 cell，并把实际生成的脚本路径和 diff/hash 打印出来。

### I12. 使用 `in dir()`、跨 cell 临时变量和静默 fallback，隐藏运行顺序依赖

- **位置**：01 系列、04、05、06、D04、D05、D06、D07、D09、D11、D13 多处；典型如 `05:cell 26/42`、`06:cell 27/29/30`。
- **影响**：用户乱序重跑时可能复用旧 kernel 中的 `results`、`_converged` 或 `_verdict_results`，得到与 fresh run 不同的输出；默认值会掩盖前置 cell 未执行。
- **改法**：PARAMS 后增加 `RUN_ID` 与状态初始化 cell；每个分析块显式检查所需变量/adata 字段，缺失即告诉用户应运行哪个 cell。不要用 `dir()` 判断科研结果是否存在；应从 `adata` 或磁盘中的版本化结果读取。

### I13. 01 per-sample 过滤影响 cell 顺序错误，报告不到过滤前后变化

- **位置**：`01_kim.ipynb:cell 28/32`、`01_nancang.ipynb:cell 28/32`、`01_yue.ipynb:cell 28/32`
- **影响**：cell 28 尝试读取 `adata.uns['qc_report_v1']`，但该对象到 cell 32 才创建；因此用户看不到承诺的每样本 before/after 过滤影响。另有乱码注释和 `pass` 占位残留。
- **改法**：在过滤前保存 per-sample counts，在过滤后立即构建 before/after/delta 表并显示、导出；删除 `dir()` 和空 `pass`。

### I14. smoke runner 会把失败流程误判为 PASS

- **位置**：`scripts/smoke_run_notebooks.py:64-181`
- **影响**：使用 `--allow-errors`；只要 executed notebook 中没有 error output 就判 PASS，不检查 nbconvert return code；没有输出的赋值 cell被判 SKIP；D01–D14 的 `output_check` 全为空；核心方法“优雅跳过”仍算 PASS。
- **改法**：移除 `--allow-errors` 或把任何 error cell/非零 return code 视为失败；按每个 notebook 声明核心产物字段和文件；区分 `executed-no-output` 与 skipped；报告每种可选方法是 success/skipped/unavailable/failed。至少用小夹具验证真实 stage contract，而不是只验证 notebook 没抛异常。

## Minor Issues

### M1. README 已明显滞后且互相矛盾

- `notebooks/README.md` 无 frontmatter，仍写不存在的 `01_*_loaded.ipynb`、错误输出名；声称所有下游彼此独立，但 D08 依赖 D04–D07，D14 消费拟时序结果。
- `07_downstream/README.md` 把 D02 写成 PyDESeq2（实际 R DESeq2）、D12 写成 LIANA（实际 CellChat/CellPhoneDB）、D14 写成 tradeSeq（实际 UnivariateSpline R2），并把 CytoTRACE2 写成 R（notebook 是 Python 包路径）。

### M2. Stage 编号与文件名/输出名混乱

- 多个 D notebook 标题仍写 07、08、10b、14、15、16；输出文件和 `uns['stage']` 使用旧编号。生物医学用户无法仅凭文件名判断依赖链。
- 建议所有界面统一使用 `D01...D14`，旧编号只放迁移说明。

### M3. 06 有明显复制合并残留

- `06_annotated.ipynb:cell 23/25` 同一条件和注释重复三次。虽不直接改变结果，但降低审查可信度，应清理。

### M4. deprecated notebook 的指引已过期

- `_deprecated/01_loaded.ipynb` 说由 `01_*_loaded.ipynb` 替代，当前实际命名不是如此；它仍展示 `read_with_manifest` 旧入口，容易被学生误用。
- 建议 deprecated 目录默认不出现在主导航，并在开头给准确替代路径，不再保留可直接 Run All 的旧代码。

### M5. 安装提示会诱导用户在错误环境直接 `pip install`

- 多个 notebook 输出 `pip install ...` 或降级 numpy 的命令，但项目要求 conda 环境隔离和精确 pin。
- 建议只指向对应 `environment-*.yml` / setup 脚本，并打印当前环境名与预期环境，不鼓励 notebook 内临时装包。

### M6. 输出命名会覆盖同参数重跑结果

- 大量图和 CSV 使用固定文件名；改参数重跑会覆盖旧结果，而 h5ad 版本号常仍为 v1。
- 建议输出目录包含 `RUN_ID` 或参数短 hash，并生成 `run_manifest.yaml` 记录输入文件 hash、参数、包版本、Git commit 和时间。

## 逐 Notebook 覆盖摘要

| Notebook | 审查摘要 |
|---|---|
| 00 | 设计期与运行期分离思路正确；存在 LLM 隐私闸门缺失和 mapping 方向矛盾 |
| 01_kim / nancang / yue | 参数集中、QC 图较丰富；SoupX 后 QC 陈旧、doublet 只标不处理、per-sample 报告顺序错误 |
| 01_nowicki | 已明确跳过重复 QC；但 normalized X 被当作可与 raw counts 合并的矩阵，且 QC 指标与其他来源不可比 |
| 01_template_10x | 教学结构较完整；应同步正式 01 的修复，并避免模板与生产 notebook 两套逻辑继续漂移 |
| 02 | 交集/平衡/metadata 诊断较好；表达尺度合并与 preprocessing 并集逻辑为核心阻断问题 |
| 03 | HVG sweep、排除和诊断较丰富；counts 契约在混合输入下失效，参数选择结果默认取“最后组合”需更显式 PI 选择 |
| 04 | 方法比较与参数说明较好；counts fallback bug、无效参数、随机性和默认多方法成本需修复 |
| 05 | 多分辨率、marker/QC/compartment 诊断方向正确；巨型指标 cell 和跨 cell `results` 隐状态需拆分 |
| 06 | 多证据界面有价值；最终标签闸门违反 PI 判断权，是阻断问题 |
| 06b | 逐簇报告符合使用场景；306 行主循环难调试，LLM 叙述需隐私/失败状态和输入证据快照 |
| 06c | subset 逻辑和标签回流解释清楚；大量硬编码参数未与 03–06 同步，final 标签仍缺显式确认闸门 |
| D01 | marker DEG 可作描述性注释；疾病 DEG 不应作正式推断 |
| D02 | 当前统计模型错误地混合 cell type，必须重写 |
| D03 | 方法步骤可读；基因位置逻辑隐藏在 src，参考标签依赖精确字符串匹配，需参考细胞 QC |
| D04 | 多指标 root 界面较丰富；stem marker raw mean 有深度混杂，NaN 标签仅 warning 不阻断 |
| D05 | R 桥可见；204 行内联脚本难调试，方法失败仍写正式 checkpoint |
| D06 | fate-mapping 主线可见；pseudotime 输出契约与序列化删除式 fallback 需收紧 |
| D07 | 尚属未验证原型，不应宣称正式打通 |
| D08 | 方法对齐/相关/Jaccard 思路合理；依赖关系应在主导航和 smoke 顺序中明确 |
| D09 | 两组 sample-level 检验相对清楚；多组模块是假统计，默认 quick-test 不应出现在正式结果路径 |
| D10 | ORA + pathway activity 双路线合理；输入基因/背景和矩阵尺度需修正 |
| D11 | 三步 SCENIC 叙事清楚；TF universe 和表达输入使结果当前不可信 |
| D12 | 多工具互补设计合理；CellChat 输入尺度和标签转换错误 |
| D13 | R 侧 hdWGCNA 流程可见；Python 聚类距离错误、巨型 R cell 难调试 |
| D14 | 可作探索性趋势可视化；不能称为 DEG/GAM 显著性分析 |
| deprecated 01/02 | 不应作为用户入口；替代路径和旧 API 说明需更新 |

## 面向生物医学专家的统一 Notebook 模板

建议所有 notebook 统一为以下 10 个区块。重点不是把逻辑抽进 `src`，而是把科研判断拆成可读、可调、可验证的 cell。

1. **本 stage 回答什么生物学问题**
   - 一句话研究问题；适用/不适用场景；本方法不能证明什么。
2. **输入/输出契约**
   - 输入路径、必需 `obs/var/layers/obsm`、表达尺度、上游 stage、正式输出字段。
3. **PARAMS：只放需要 PI/学生调整的参数**
   - 每个参数给单位、允许值、推荐范围、增大/减小的后果、何时不该调整。
4. **运行前检查**
   - 数据维度、样本数、每组 donor 数、缺失值、counts 是否非负近整数、标签覆盖率、内存估算、环境/设备。
5. **输入诊断图与 PI 决策点**
   - 在运行方法前先画决定参数所需的图；提供“看什么、什么情况回上游”。
6. **方法执行：每个方法独立 20–60 行 cell**
   - 可选方法用 `if METHOD_ENABLED:` 包完整区块；不要把多个方法塞进同一巨型 cell。
7. **中间结果 QC**
   - 明确成功标准、失败标准、数值范围、细胞/样本覆盖率、负对照/已知 marker、收敛图。
8. **PI 明确选择/确认 cell**
   - 自动建议只写 suggested 字段；final 字段必须有 `PI_CONFIRMED=True` 和人工文件输入。
9. **追溯与保存**
   - `RUN_ID`、输入 hash、Git commit、参数、随机种子、包版本、核心产物列表、四态状态。
10. **Stage Verdict 与下一步**
   - checkbox 不是装饰；用 assert/硬闸门保证必需项完成，列出应重跑的下游 notebook。

### PARAMS cell 示例规范

```python
# === PARAMS：PI/学生只修改本 cell ===
UPSTREAM_PATH = "results/..."
OUTPUT_DIR = "results/runs/<RUN_ID>/..."

# 科研选择：按 donor 做 pseudobulk；不能改成 cell_id
SAMPLE_COL = "donor_id"

# 最少独立生物学重复；低于 3 只允许描述性展示，不做显著性推断
MIN_DONORS_PER_GROUP = 3

# 默认关闭；只有查看完上方诊断并完成 CSV 决策后改 True
PI_CONFIRMED = False

RANDOM_SEED = 42
```

### 每个方法结尾的四态验收

```python
method_status = {
    "state": "success",  # success/skipped_by_user/unavailable/failed
    "reason": "",
    "required_outputs": ["obs[...]", "tables/...csv"],
}
assert method_status["state"] == "success", method_status
```

## 优先级路线图

### P0：先修科研正确性，暂停正式分析

1. 重构 01–03 表达尺度契约，禁止 raw/normalized 混合；为每来源记录真实 counts 可用性。
2. 重写 D02 为“逐 cell type、sample/donor 级 DESeq2”。
3. 给 06 增加 PI 显式确认闸门，禁止自动生成 final 标签。
4. 修正 D11 TF universe、D12 CellChat 输入、D14 方法命名/统计模型。
5. 给 00 增加临床隐私脱敏与外部 LLM 双确认。

### P1：把 notebook 变成真正可调、可调试的专家工作台

1. 将所有 >100 行巨型 cell 按准备/执行/解析/QC/保存拆分。
2. 清除 `in dir()` 和旧 kernel 隐状态；每块使用显式输入契约。
3. 核心分析逻辑从 `src/io.py` 回到 notebook，展开基因映射与位置注入证据。
4. 所有可选方法采用四态结果，失败不写正式 checkpoint。
5. 为每个 stage 增加参数生效摘要、方法成功标准和 PI decision cell。

### P2：建立真实的 notebook 验收体系

1. 修复 smoke runner：非零返回码/error cell/缺核心字段均失败。
2. 为每个 notebook 建小型但生物学合理的 fixture 和 contract test。
3. 对 scVI、Monocle3、CellRank2、CytoTRACE2、CellChat、hdWGCNA 在支持环境各做一次真实运行并保存基准报告。
4. 增加 fresh-kernel 顺序执行测试和关键 cell 单步重跑测试。

### P3：统一文档与追溯

1. 重写两级 README，与实际工具、路径、依赖关系和输出字段一致。
2. 统一 D01–D14 编号；deprecated 入口退出主导航。
3. 所有输出用 `RUN_ID/参数 hash`，生成 `run_manifest.yaml`，避免覆盖并支持论文复现。

## 最终判断

当前 notebook 的“解释性外壳”已经比普通科研脚本友好，但**科研执行契约还没有达到可放心交给生物医学专家的程度**。下一轮优化不应优先继续增加分析方法，而应先把 P0 的数据尺度、统计单位、人工确认和方法输入四类问题修正；否则越丰富的下游模块只会放大上游错误。
