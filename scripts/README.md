---
title: 脚本目录说明
tags: [scripts, guide]
created: 2026-07-10
updated: 2026-07-12
---

# scripts — 工具脚本

本目录存放流水线主线（`notebooks/`）之外的辅助脚本,包括环境搭建、测试数据
生成、以及端到端冒烟运行器。

## 脚本清单

| 脚本 | 语言 | 用途 |
|------|------|------|
| `make_test_subset.py` | Python | 从完整数据集抽样生成小型 h5ad,供 pytest 集成测试和 CI 使用 |
| `smoke_run_notebooks.py` | Python | 端到端冒烟运行器：依次 nbconvert --execute 主线（01 到 06）与全部下游（D01 到 D14）notebook,逐 cell 记录 pass/fail 并产出 JSON 汇总,验证 h5ad 输入输出契约不断裂。不是 pytest 用例（不被 `tests/` 收集）,命令行直接运行 |
| `env_parity.py` | Python | 跨平台环境一致性诊断：对比两台机器的 conda 环境快照并列出差异（ADR-0010）|
| `trace_downstream.py` | Python | 扫描 07_downstream 各 notebook,提取其读取与写入的 obs 列,追踪列依赖关系 |
| `deseq2_contrast.R` | R | `D02_pseudobulk_deg.ipynb` 调用的 DESeq2 差异表达对比 |
| `soupx_run.R` | R | `01` QC 阶段可选的 SoupX 环境 RNA 污染校正 |
| `install_r_github_pkgs.R` | R | 安装需从 GitHub 源码编译的 R 包（Monocle3、CytoTRACE2、hdWGCNA）,这些包不在 conda 渠道 |
| `setup_cuda.sh` | Bash | 在有 NVIDIA CUDA 的 Linux 机器上,把 CPU 版 PyTorch 替换为 CUDA 版（ADR-0013）|
| `prefetch-gitleaks.sh` | Bash | 预取 gitleaks 检测工具,加速 pre-commit 首次运行 |
