---
title: "scRNA-seq整合分析框架"
id: "scrna-integration-framework"
type: "research"
status: active
phase: analysis
priority: "high"
tags: [single-cell, scRNA-seq, integration, framework, bioinformatics, multi-disease]
created: "2026-06-05"
updated: "2026-07-10"
external_path: ""
visibility: public
repo:
  url: git@github.com:kaisermoon/scrna-integration-framework.git
  default_branch: main
  protected_branches: [main]
  ci_required_checks: [test, lint]
  pr_size_limit: {files: 15, lines: 3000}  # 认知复杂度用「改动 cell 数 ≤ 30」控制；notebook JSON 膨胀不计入，详见 CLAUDE.md
  reviewer_profile: default
---

# scRNA-seq整合分析框架

## 项目概述

构建一套**跨病种、多源 scRNA-seq 数据整合分析框架**：支持原始 cellranger 矩阵 / h5ad / RData 等多源接入，标准化预处理与质控，多方法去批次（Harmony / scVI / scANVI 等）与多方法注释（marker / LLM / scANVI / AUCell / UCell）并存比较，下游模块化扩展（拟时序、GRN 等）。每步骤模块化、可替换、可迭代回跑，以 Jupyter notebook 作为主要交付形式，支持研究者 - AI 协作不断打磨参数。

## 核心设计原则

1. **多源接入**：cellranger / h5ad / RData / 任意 obs schema → 统一为 CellxGene-compatible 格式
2. **模块化与可替换**：每步骤约束输入输出（h5ad），任意环节可换方法做对比
3. **迭代回跑**：注释发现问题 → 回退到 HVG / embedding / clustering 重跑，机制内建
4. **多方法并存**：去批次、注释、评分等关键决策步骤允许多方法同时跑，结果存于 obs 列做交叉比对
5. **内存自律**：稀疏矩阵 + 避免对象复制，面向海量整合
6. **学术追新**：新方法发表 → 快速插入流水线，不重写下游

## 相关项目

- [[钟子劭个人工作全景]]
- 上游知识库：[[scRNA-seq方法学]]（待建）

## 外部目录

- **路径**：项目目录本身就是 GitHub 仓库根，无外部路径
- **大数据存放**：原始数据放 `data/raw/`（gitignore）；分析对象放 `results/data/`（gitignore）
- **早期框架参考**：`/Users/zhongzishao/Works/GCPL_scRNA/`（已部分复制到 `references/legacy-GCPL/`）

## 结构现状

代码组织遵循「src 与 notebook 边界铁律」（见项目级 `CLAUDE.md`）：`src/scrna_integration/` 只收无科研价值的技术管道，一切分析逻辑留 notebook cell。

### src 技术管道模块

| 模块 | 文件 | 职责 |
|------|------|------|
| 启动脚手架 | `src/scrna_integration/bootstrap.py` | 定位项目根 + BLAS 线程设置，notebook 顶部一行调用 |
| 平台检测收口 | `src/scrna_integration/platform.py` | OS / CUDA / Rscript 路径检测，`detect_device()` 设备选择单点收口（ADR-0013） |
| LLM 配置与调用 | `src/scrna_integration/llm_config.py` | 06 系列注释的多模型直连配置与标准调用 |
| 多格式读取 | `src/scrna_integration/io.py` | 基因 ID 双向同步 + 基因组位置注入 + batch 键诊断（纯技术管道） |

### notebooks 分析流水线

| 阶段 | notebook | 说明 |
|------|----------|------|
| 前置 | `notebooks/00_propose_obs_manifest.ipynb` | obs 映射清单辅助生成 |
| 每来源接入 + QC | `notebooks/01_per_dataset/` 四个格式模板（`01_template_10x_mtx` / `01_template_10x_h5` / `01_template_h5ad` / `01_template_counts_matrix`） | 按**输入数据格式**命名（数据集属下游项目资产，格式才是框架骨架）；接新数据集复制对应格式模板改 PARAMS；字段映射/QC/断言在 cell |
| 合并 | `notebooks/02_merged.ipynb` | concat + 合并质控 |
| 归一化 / HVG | `notebooks/03_normalized.ipynb` | 归一化状态判定 + HVG |
| 降维去批次 | `notebooks/04_embedded.ipynb` | PCA/Harmony/scVI/scANVI 多方法并存比较 |
| 聚类 | `notebooks/05_clustered.ipynb` | 多分辨率 Leiden + 多指标推荐 |
| 注释 | `notebooks/06_annotated.ipynb` + `06b_per_cluster` + `06c_subset` | 多方法交叉比对 + LLM 判决 |
| 下游 | `notebooks/07_downstream/`（`D01`–`D14`） | 各模块以 `06_annotated.h5ad` 为输入、彼此无执行先后、可单独运行 |

## 数据集（已收集，部分）

参见 `/Users/zhongzishao/Works/GCPL_scRNA/data/`：Kim 2025 / Nancang 2025 / Nowicki-Osuch 2023 / Tsubosaka 2023 / Yue 2025 等多个 GEO 数据集，整合方向为胃癌前病变 / 胃癌相关单细胞数据。后续将正式纳入框架的首个验证案例。
