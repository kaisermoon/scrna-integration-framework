---
title: 脚本目录说明
tags: [scripts, guide]
created: 2026-07-10
updated: 2026-07-10
---

# scripts — 工具脚本

本目录存放流水线主线（`notebooks/`）之外的辅助脚本,包括环境搭建、测试数据
生成、以及部分下游分析的无界面运行器。

## 命名约定

带下划线前缀的脚本（`_run_*.py`）是某个 notebook 核心逻辑的无界面（headless）
运行器。它们把 notebook 里的计算步骤抽出来单独成脚本,用途是：scVI 训练、
大规模 R 分析等步骤通过 nbconvert 端到端运行会超时,改用这些脚本在命令行直接
运行以完成验证。它们与对应 notebook 的逻辑保持一致,不是独立功能。

## 脚本清单

| 脚本 | 语言 | 用途 |
|------|------|------|
| `make_test_subset.py` | Python | 从完整数据集抽样生成小型 h5ad,供 pytest 集成测试和 CI 使用 |
| `run_pipeline_test.py` | Python | 端到端冒烟测试：按顺序运行 01 到 06 核心 stage,验证 h5ad 输入输出契约不断裂 |
| `env_parity.py` | Python | 跨平台环境一致性诊断：对比两台机器的 conda 环境快照并列出差异（ADR-0010）|
| `trace_downstream.py` | Python | 扫描 07_downstream 各 notebook,提取其读取与写入的 obs 列,追踪列依赖关系 |
| `_run_07_deg.py` | Python | `07_deg.ipynb` 的无界面运行器 |
| `_run_09_cnv.py` | Python | `09_cnv.ipynb`（infercnvpy）的无界面运行器 |
| `_run_nb_16_trajectory_de.py` | Python | `16_trajectory_de.ipynb`（tradeSeq）的无界面运行器 |
| `deseq2_contrast.R` | R | `08_pseudobulk_deg.ipynb` 调用的 DESeq2 差异表达对比 |
| `soupx_run.R` | R | `01` QC 阶段可选的 SoupX 环境 RNA 污染校正 |
| `install_r_github_pkgs.R` | R | 安装需从 GitHub 源码编译的 R 包（Monocle3、CytoTRACE2、hdWGCNA）,这些包不在 conda 渠道 |
| `setup_cuda.sh` | Bash | 在有 NVIDIA CUDA 的 Linux 机器上,把 CPU 版 PyTorch 替换为 CUDA 版（ADR-0013）|
| `prefetch-gitleaks.sh` | Bash | 预取 gitleaks 检测工具,加速 pre-commit 首次运行 |
