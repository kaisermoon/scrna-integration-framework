# ADR-0011: Per-Dataset QC Notebooks over Monolithic Load-then-QC

- **Status**: Accepted
- **Date**: 2026-06-11
- **Supersedes**: None (extends existing stage structure)

## Context

原 01_loaded 读单个 manifest，02_qcd 在合并后的 AnnData 上循环做 per-dataset QC。这导致：

1. 无法为每个数据集独立调 QC 阈值（组织活检 vs 类器官的 MT% 基线完全不同）
2. 研究者无法独立查看每个数据集的质量分布后再决定合并
3. 合并时机不明确（SPEC 中未定义 `read_with_manifest` 多次调用如何 concat）
4. QC 参数空间过小（仅 4 个固定阈值），缺乏自适应策略

## Decision

**拆分 01-02 为 per-dataset notebooks + 显式 merge notebook**：

```
notebooks/01_per_dataset/{dataset}.ipynb  →  per-dataset h5ad
notebooks/02_merged.ipynb                 →  合并 h5ad（03 的 upstream）
```

每个 per-dataset notebook 独立完成：`read_with_manifest` → 完整 QC（自适应阈值、双细胞、SoupX、细胞周期等）→ checkpoint。

02_merged 做显式 `anndata.concat` + 跨数据集诊断。

### 关键设计选择

1. **MAD-based 自适应阈值**替代固定阈值（N_MAD 倍数是真旋钮）
2. **Notebook IS the QC report**——无独立报告层，所有诊断图在 notebook 内即时输出
3. **标记不过滤**原则——血红蛋白/应激/doublet 标记为 obs 列，不自动移除，PI 在查看后决定
4. **Scalar-or-Sweep 双模参数**（03-05）——单值直接跑，列表自动 sweep + 对比图

## Consequences

- 新增 4 个 per-dataset notebook + 1 个 merge notebook
- 旧 01_loaded + 02_qcd 移入 `notebooks/_deprecated/`
- `read_with_manifest()` 函数不变（它本就是 per-manifest 设计）
- 下游 03-15 合约不变（02_merged 输出格式等同旧 02_qcd）
- 每新增一个数据集 = 新增一个 per-dataset notebook + 修改 02_merged 的 PER_DATASET_PATHS

## Alternatives Considered

1. **单 notebook + manifest 列表循环**：调参不便，无法独立看各数据集
2. **配置文件驱动的 per-dataset QC 参数**：违反 ADR-0001 薄框架，参数应在 PARAMS 块可见
3. **Snakemake 编排**：过重，不符合 notebook 交互工作台定位
