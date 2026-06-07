---
title: "执行计划：scRNA-seq整合分析框架"
updated: "2026-06-07"
phase: planning
---

## 当前方向

第三轮 grilling（2026-06-07，对照原始构思 + references 真实代码 GCPL/student-code）触发 6 处修订，落进 ADR-0006~0008 + SPEC/CONTEXT 同步：

1. manifest 用 YAML，确立"写一次数据事实→manifest / 要调旋钮→PARAMS"边界（ADR-0006）
2. R 桥接按工具分流（纯 Python / rpy2 轻量 / subprocess 重型），推翻 ADR-0002（ADR-0007）；InferCNV 改回纯 Python infercnvpy
3. stage6 注释 4 默认同跑 + CellTypist 候选（多方法=交叉验证）
4. stage5 纯 Leiden+sweep，聚类方法走 obs 列并存扩展模式（不要 ACDC 默认）
5. stage4 所有 embedding 平级，UMAP 目测为主决策 + 指标辅助
6. student-code 下游全部吸收进规划但按规范重写（ADR-0008），stage7 扩展每模块独立 PR

PR 计划重排：老 PR-3 拆为 PR-3a/3b/3c（主线到注释）+ PR-4（第一波结果）+ PR-5（marker）+ PR-6~11（stage7 扩展）。下一步：进入实施期，委派 coder 起 PR-0a。

## 已确定的架构层（详见 SPEC.md + 5 条 ADR；术语见 CONTEXT.md）

1. **定位**：Stage 1 个人科研基础设施（PI + 学生 + AI agent），逐步可升级 software paper / PyPI
2. **薄框架到底**：只 2 个 Python 函数 `read_with_manifest` + `sweep`；其他全部回归 scanpy / anndata / pandas / yaml 原生
3. **stage 切片 = 文件命名约定**：1 / 2 / 3 / 4 / 5 / 6 / 6.5 / 7（含 9 个 downstream 模块），不是 Stage 类，没有 entry function
4. **多方法并存**：用 anndata 自然 `obsm` / `obs` 多 slot，scanpy 默认命名（`X_pca` / `leiden_res_0.5` 等）
5. **回跑机制 = 文件名版本号 + `adata.uns["status"]` PI 直接写**
6. **rpy2 + anndata2ri** 直接用，不包装 `_r_bridge`，notebook 里像 legacy-GCPL 那样写
7. **Sweep harness**：`sweep(fn, adata, candidates, scorer, output_dir) → DataFrame`，scorer 是普通函数
8. **obs Schema**：CellxGene 软对齐，三层约定写在文档；`read_with_manifest` 仅做 manifest schema 校验，不做 obs 校验
9. **Manifest** schema 完整保留（含 obs_mapping / value_mapping / clinical_metadata / ontology / project_specific / original_annotations / preprocessing_done / qc_overrides）
10. **Marker 库**：`references/markers/*.csv`，pandas 直接读
11. **Disease ontology**：`references/disease_ontology/*.yaml`，yaml.safe_load 直接读，PI 自己 walk dict
12. **QC heterogeneity**：stage 2 notebook (`stage2_qcd.ipynb`) 内 if/else 处理 manifest 的 `preprocessing_done` + `qc_overrides`
13. **内存纪律**：CSR / inplace / del / lzf / float32——5 条约定写文档，code review 卡
14. **Stage 6 注释 cross-method comparison**：`stage6_annotated.ipynb` 内的 5 方法 + LLM verdict + PI 拍板，**没有 si.report.stage6_annotation**
15. **Per-cluster 深度报告**：`stage6_per_cluster.ipynb` + 普通 for-loop，**没有 plugin 系统**

## 任务表 — 渐进推进 PR 计划

> 第三轮 grilling 重排（2026-06-07）：老 PR-3 的"12 notebook 一次端到端"过重，拆为主线（到注释）→ 第一波结果 → marker 库 → stage7 扩展（每模块独立 PR，吸收 student-code 技术点并按本项目规范重写，见 ADR-0008）。

### 第一阶段：框架骨架（PR-0a → PR-2）

| # | PR | 范围 | 大小 | 状态 | 验收 |
|---|------|------|------|------|------|
| 1 | PR-0a | pyproject 依赖（scanpy / anndata / scvi-tools / harmonypy / scrublet / **infercnvpy** / **cellrank**（CytoTRACE）/ **decoupler**（pseudobulk）/ rpy2 / anndata2ri / openrouter SDK / mLLMCelltype / mygene / **sccoda** 等）+ environment.yml + environment-r.yml（含 SoupX / DESeq2 / Monocle3 / UCell / hdWGCNA 等 R 包，按 stage 分组可选装） | ≤10/400 | pending | smoke test 通过 + 双环境可装 |
| 2 | PR-1 | `src/scrna_integration/io.py`：`read_with_manifest()` 完整实现（10x_mtx + h5ad + RDS via rpy2 + obs_mapping + value_mapping + clinical_metadata join + original_annotations 重命名 + 基因 ID 双向同步 + species 校验 + raw matrix 路径记录 + preprocessing_done/qc_overrides schema 校验）+ 5 个 GCPL manifests（最小必填优先，仅异构数据集写可选块）+ 最小 gastric ontology（5-10 节点） | 申请 size 例外 | blocked PR-0a | **GCPL 5 数据集 stage 1 端到端跑通** + 每个 dataset 出一个 stage1_loaded_v1.h5ad，含 `var["ensembl_id"]` + `adata.uns["species"]` + `adata.uns["raw_matrix_path"]`（如有） |
| 3 | PR-2 | `src/scrna_integration/sweep.py`：`sweep()` 完整实现 + `src/scrna_integration/scorers.py`：常用 scorer（QC 平衡、`integration_metrics` scIB suite、silhouette + ARI clustering_metrics、annotation concordance）+ `src/scrna_integration/markers.py`：`load_markers()`（按 ADR-0005）+ `__init__.py` re-export 三函数 | ≤10/400 | blocked PR-1 | unit test：sweep wrap 任意 scanpy/scvi 函数；load_markers 三种 role 模式；integration_metrics 在合成数据上输出合理 |

### 第二阶段：主线 notebook 跑到注释（PR-3a → PR-3c）

| # | PR | 范围 | 大小 | 状态 | 验收 |
|---|------|------|------|------|------|
| 4 | PR-3a | 主线前段 notebook：`stage1_loaded` / `stage2_qcd`（manifest 驱动 QC skip + scrublet + SoupX rpy2 按需读 raw matrix）/ `stage3_normalized`。每个含 PARAMS cell + scanpy 原生 + 内存纪律 idiom + self-check + 结尾 del/gc | ≤10/400 | blocked PR-2 | **GCPL 5 数据集端到端跑通 stage 1→3**，出 stage3_normalized_v1.h5ad |
| 5 | PR-3b | 主线中段 notebook：`stage4_embedded`（多 embedding 平级 + UMAP 三上色目测 + sweep integration_metrics）/ `stage5_clustered`（多分辨率 Leiden + sweep + 聚类扩展槽）| ≤10/400 | blocked PR-3a | **跑通 stage 4→5**，≥1 次 stage4 重跑触发 v2 + `adata.uns["status"]="promoted"`；UMAP 多 embedding 对比图产出 |
| 6 | PR-3c | 注释段 notebook：`stage6_annotated`（4 默认注释 + scANVI 有图谱时 + CellTypist 注释 cell + cross-method comparison + LLM verdict）/ `stage6_per_cluster` / `stage6_5_subset` | 申请 size 例外 | blocked PR-3b + PI revoke/新建 OpenRouter key | **跑通 stage 6 + 6.5**，出 cell_type_final_v1 + 每簇 LLM verdict markdown |

### 第三阶段：GCPL 第一波生物学结果（PR-4）

| # | PR | 范围 | 大小 | 状态 | 验收 |
|---|------|------|------|------|------|
| 7 | PR-4 | stage 7 核心 3 模块：`stage7/deg.ipynb`（rank_genes_groups）/ `stage7/pseudobulk_deg.ipynb`（decoupler + DESeq2 subprocess）/ `stage7/cnv.ipynb`（infercnvpy 纯 Python，参考 student `4.2`/`4.3` 重写） | ≤10/400 | blocked PR-3c | **出至少一项可看生物学结果**：CAG→IM 跨阶段 DEG list / 区分肿瘤 vs 正常 CNV 图 |

### 第四阶段：marker 库 + 文档（PR-5）

| # | PR | 范围 | 大小 | 状态 | 验收 |
|---|------|------|------|------|------|
| 8 | PR-5 | `references/markers/` 初始库（gastric_epithelial / immune / stromal）+ marker 使用约定（gene 存在性检查 idiom）+ 框架 README 更新 + ADR 索引页 | ≤10/400 | blocked PR-4 | PI review 通过 |

### 第五阶段：stage 7 扩展模块（PR-6+，每模块独立 PR，吸收 student 技术并重写）

| # | PR | 范围 | 大小 | 状态 | 验收 |
|---|------|------|------|------|------|
| 9 | PR-6 | `stage7/pseudotime.ipynb`：转录组熵（numpy）+ CytoTRACE（cellrank）+ 多指标 root 识别（Z-score 综合，参考 `4.4`/`4.5`）+ Monocle3 轨迹（subprocess，参考 `4.3`/`11.2`）| 申请 size 例外 | blocked PR-5 | GCPL 上皮谱系跑通，出拟时序 + 定根结果 |
| 10 | PR-7 | `stage7/abundance.ipynb`：scCODA 贝叶斯组成差异 + Mann-Whitney/Cliff's delta/效应量（参考 `11_all_*`）| ≤10/400 | blocked PR-5 | 跨疾病阶段细胞组成差异结果 |
| 11 | PR-8 | `stage7/pathway.ipynb`：GSEApy / decoupler / Reactome | ≤10/400 | blocked PR-5 | 通路富集结果 |
| 12 | PR-9 | `stage7/grn.ipynb`：pySCENIC | 申请 size 例外 | blocked PR-5 | TF 调控网络结果 |
| 13 | PR-10 | `stage7/cell_communication.ipynb`：CellChat / CellPhoneDB | ≤10/400 | blocked PR-5 | 配体-受体通讯结果 |
| 14 | PR-11 | `stage7/gene_modules.ipynb`：hdWGCNA（subprocess）| ≤10/400 | blocked PR-5 | 共表达模块结果 |

**状态值**：`pending` | `in_progress` | `done` | `failed` | `blocked`

PR-0a → PR-1 → PR-2 → PR-3a → PR-3b → PR-3c → PR-4 → PR-5 严格串行（主线 + 第一波结果）；PR-6+ 都 base 在 PR-5 之后，彼此独立，按 GCPL 实际分析需求排序，不强制全做。

## 阻塞

- PR-1 size 远超 ≤10/400 上限（read_with_manifest 涵盖 IO/schema/clinical join/RDS bridge 全栈 + 5 manifests + 1 ontology），需独立 PR 申请 size 例外
- PR-3c / PR-6 / PR-9 size 可能超限（注释段三 notebook / 拟时序四技术 / SCENIC），需 size 例外
- 5 个 GCPL 数据集中 Yue_SSK_2025 / Tsubosaka_2023 等的具体格式需在 PR-1 实施时确认；如未在 `data/` 物理就位，需先从 `~/Works/GCPL_scRNA/data/` 复制或软链入项目 `data/`
- PR-3c stage 6 notebook 实施前需 PI revoke 旧 OpenRouter key 并新建一个，存放 `.env`
- stage 7 扩展模块（PR-6+）吸收 student-code 时必须按 ADR-0008 重写，禁止整段复制；reviewer 卡 Windows 硬编码路径 / `!pip install` cell / 800 行单体脚本 / 复制粘贴模板族

## 关键决策

| 日期 | 决策 | 理由 |
|------|------|------|
| 2026-06-05 | 项目工程形态：`src/scrna_integration/` + `notebooks/` + `tests/` | 既要可发布、也要可复现 |
| 2026-06-05 | 框架以 anndata/scanpy 为主干，CellxGene schema 为 obs 标准 | 与社区主流对齐 |
| 2026-06-05 | 80% 直接调 scanpy（ADR-0001） | scanpy API 是事实标准 |
| 2026-06-05 | rpy2 + anndata2ri 全栈统一（ADR-0002） | rpy2 在 QC stage 2 SoupX 已无法回避 |
| 2026-06-05 | Sweep 永不 auto-pick | 与 SOUL.md "核心判断权不外包" 对齐 |
| 2026-06-05 | 注释 cross-method comparison 砍掉自动多数投票 | 多 LLM 共识可能共享系统性偏差 |
| 2026-06-05 | 验证策略：直接 GCPL 5 数据集 mini 子集 | 真实场景比 sanity 更早暴露问题 |
| 2026-06-05 | PR 拆分按 stage 流程顺序 | 框架质量持续端到端验证 |
| 2026-06-06 | 朴素优先于 plugin 系统（ADR-0003） | 第一轮 grilling 反复出现"加 decorator + registry"的过度工程倾向 |
| 2026-06-06 | 框架大瘦身：从 11 个 si.* API 砍到 2 个函数（ADR-0004） | 第二轮 grilling 后 PI 反对所有 wrapper；legacy GCPL 全是 scanpy 原生 |
| 2026-06-06 | `load_markers` 加为合法第 3 函数（ADR-0005） | PI 实战经验确认 boilerplate 真重复 + role 语义集中 |
| 2026-06-06 | Notebook 不要 `_template_` 前缀，直接可跑代码 | PI 第三轮反馈：模板增加 copy/rename 仪式无价值 |
| 2026-06-06 | 项目目的=生物学发现，不是方法学论文 | 节 1 PI 重申；删除 Stage 2 software paper 期刊清单 + Stage 3 PyPI 目标 |
| 2026-06-06 | 基因 ID 双向同步必做（symbol + ensembl 同时存在） | 多源整合的硬需求；GCPL 5 数据集就有 symbol/ensembl 混用 |
| 2026-06-06 | 物种声明必填，目前只接受 human | 框架不预留跨物种；如未来加 mouse 模型走新 ADR |
| 2026-06-06 | stage 6.5 subset analysis 作为正式阶段 | 学生代码已用此模式；PI 实际研究中 T 细胞 / SPEM 亚群分析常需 |
| 2026-06-06 | stage 4 integration QC 嵌 sweep 报告，不另立 stage 4.5 | 与 sweep 复合 scorer 模式对齐 |
| 2026-06-06 | stage 7 拆 9 模块；首期 PR-3 只做核心 3 模块（DEG / pseudobulk DEG / CNV）；其他 6 模块按 GCPL 需要 PR-5+ 逐一加 | 节 1 项目目的=出生物学发现，首期 milestone 应对齐 GCPL 第一波结果 |
| 2026-06-07 | manifest 用 YAML，区分"写一次数据事实→manifest / 要调旋钮→PARAMS"（ADR-0006） | 第三轮 grilling：对照原始构思"YAML 难懂"的表面矛盾，PI 厘清数据配置配好不动、参数才反复调 |
| 2026-06-07 | R 桥接按工具分流：纯 Python（infercnvpy/CytoTRACE）/ rpy2（SoupX 轻量）/ subprocess Rscript（Monocle3/UCell/DESeq2 重型），推翻 ADR-0002 全栈 rpy2（ADR-0007） | student-code 真实下游全是 subprocess；rpy2 重对象转换脆、抗升级差；InferCNV 其实纯 Python 被 SPEC 错标 |
| 2026-06-07 | stage6 注释 4 默认同跑（marker/LLM/基因集/scANVI）+ CellTypist 候选；多方法为交叉验证非冗余 | PI："只有多几个相互比对才能做出最准确注释"；scANVI 要、CellTypist 未必 |
| 2026-06-07 | stage5 默认纯 Leiden + sweep；ACDC 等其他聚类走"写 obs 列并存"扩展模式，不进默认依赖 | PI 未跑通 ACDC（太慢），但需聚类方法可替换接口；Leiden 有时定不准群数 |
| 2026-06-07 | stage4 所有 embedding 平级，无主力；决策=UMAP 三上色目测（主）+ integration_metrics；census 预训练 scVI/scANVI 均可选 | PI："没有方法在每个数据集都好"，全跑+拉 UMAP 目测+指标再挑；未来还要加 scCARFT 等 |
| 2026-06-07 | student-code 下游全部吸收进规划，但按本项目规范重写、禁整段复制（ADR-0008）；stage7 扩展每模块独立 PR | PI：要做这些且参考其代码，但规范重写；分阶段不要一次太重 |
| 2026-06-07 | 老 PR-3（12 notebook 一次端到端）拆为 PR-3a/3b/3c（主线到注释）+ PR-4（第一波结果）+ PR-5（marker）+ PR-6~11（stage7 扩展） | PI：分不同阶段，不要一次太重 |

## 最近 insight

- 框架价值在"多源 IO + manifest + 临床 metadata join"和"sweep 循环 + 报告拼装"两处真空白；其他全部回归 scanpy/anndata/pandas/yaml 原生
- LLM judge 在 stage 6 注释 cross-method comparison + sweep recommendations 仍是核心能力，但实现位置改为 stage notebook 内的 OpenRouter 直接调用，不包装为 framework 函数
- Notebook 模板从"次要交付物"升级为"框架的主要标准化单位"——PI 通过模板传播分析流程，而非通过 API
- 反思：第一轮 grilling 我（主 Agent）反复推过度工程方案（包装 scanpy / register_run / decorator plugin / lineage system），PI 多次纠偏。这是 LLM agent 的训练分布偏差——对"完整框架"的本能倾向。code_reviewer 应在 PR 阶段持续卡这条
