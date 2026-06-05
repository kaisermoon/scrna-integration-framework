---
title: "scRNA-seq整合分析框架"
id: "scrna-integration-framework"
type: "research"
status: active
phase: planning
priority: "high"
tags: [single-cell, scRNA-seq, integration, framework, bioinformatics, multi-disease]
created: "2026-06-05"
updated: "2026-06-05"
external_path: ""
visibility: public
repo:
  url: git@github.com:kaisermoon/scrna-integration-framework.git
  default_branch: main
  protected_branches: [main]
  ci_required_checks: [test, lint]
  pr_size_limit: {files: 10, lines: 400}
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

## 子模块（计划中）

| 模块 | 目录 | 状态 |
|------|------|------|
| 数据接入与标准化 | `src/io/` + `notebooks/01_*` | 待启动 |
| 预处理与质控 | `src/qc/` + `notebooks/02_*` | 待启动 |
| 标准化 / HVG | `src/preprocessing/` + `notebooks/03_*` | 待启动 |
| 降维与去批次（多方法并存） | `src/embedding/` + `notebooks/04_*` | 待启动 |
| 聚类（参数扫描） | `src/clustering/` + `notebooks/05_*` | 待启动 |
| 注释（多方法交叉） | `src/annotation/` + `notebooks/06_*` | 待启动 |
| 差异基因 / 通路 / 基因集评分 | `src/downstream/dge/` + `notebooks/07_*` | 待启动 |
| 拟时序 | `src/downstream/pseudotime/` | 待启动 |
| GRN | `src/downstream/grn/` | 待启动 |

## 数据集（已收集，部分）

参见 `/Users/zhongzishao/Works/GCPL_scRNA/data/`：Kim 2025 / Nancang 2025 / Nowicki-Osuch 2023 / Tsubosaka 2023 / Yue 2025 等多个 GEO 数据集，整合方向为胃癌前病变 / 胃癌相关单细胞数据。后续将正式纳入框架的首个验证案例。
