"""per_dataset notebook 输出字段权威约束。

每个 01_per_dataset notebook 产出的 h5ad 必须满足本模块定义的列名、类型与取值约束，
下游 02_merged 及测试依赖此契约。

约束性质：纯机械校验，不含科研判断（如阈值合理性、双细胞率判定）。
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

# ── obs 列 ───────────────────────────────────────────────────────────────
# 双细胞三态列（所有来源必须产出，Kim/Nancang 标准）
DOUBLET_CLASS_COL = "doublet_class"          # Category: "singlet"/"uncertain"/"doublet"
DOUBLET_INCLUDE_COL = "doublet_include"       # bool: 是否纳入下游
DOUBLET_PREDICTED_COL = "predicted_doublet"  # bool: 高置信 doublet（scrublet 原始）
DOUBLET_SCORE_COL = "doublet_score"          # float: scrublet 分数

# QC 指标列（经 scanpy.pp.calculate_qc_metrics 产生）
QC_N_GENES_COL = "n_genes_by_counts"
QC_TOTAL_COUNTS_COL = "total_counts"
QC_PCT_MT_COL = "pct_counts_mt"

# ── uns 键 ────────────────────────────────────────────────────────────────
DOUBLET_CONTRACT_KEY = "doublet_contract"    # dict: scrublet 参数与结果摘要
EXPRESSION_CONTRACT_KEY = "expression_contract"  # dict: counts 来源契约

# ── layers ────────────────────────────────────────────────────────────────
LAYER_COUNTS = "counts"                      # 原始整数 counts（只读，不得覆盖）
LAYER_COUNTS_SOUPX = "counts_soupx"         # SoupX 校正后 counts（若运行 SoupX）

# ── 公开导出 ──────────────────────────────────────────────────────────────
__all__ = [
    "DOUBLET_CLASS_COL",
    "DOUBLET_INCLUDE_COL",
    "DOUBLET_PREDICTED_COL",
    "DOUBLET_SCORE_COL",
    "QC_N_GENES_COL",
    "QC_TOTAL_COUNTS_COL",
    "QC_PCT_MT_COL",
    "DOUBLET_CONTRACT_KEY",
    "EXPRESSION_CONTRACT_KEY",
    "LAYER_COUNTS",
    "LAYER_COUNTS_SOUPX",
    "validate_per_dataset_output",
]

# ── 验证函数 ──────────────────────────────────────────────────────────────

_VALID_DOUBLET_CLASSES = frozenset({"singlet", "uncertain", "doublet"})


def validate_per_dataset_output(adata: Any) -> dict[str, Any]:
    """检查 per_dataset notebook 产出的 h5ad 是否符合下游契约。

    纯机械校验，不做科研判断。返回报告字典，函数本身不 raise，
    由调用者根据 ``passed`` 决定行为。

    Parameters
    ----------
    adata : anndata.AnnData
        待校验的 AnnData 对象。

    Returns
    -------
    dict
        {"passed": bool, "errors": [str], "warnings": [str]}
    """
    errors: list[str] = []
    warnings: list[str] = []

    # 1. layers["counts"] 存在 + 稀疏矩阵 + 非负 + 近整数
    if LAYER_COUNTS not in adata.layers:
        errors.append(f"layers['{LAYER_COUNTS}'] 缺失：per_dataset 必须建立原始 counts layer")
    else:
        _counts = adata.layers[LAYER_COUNTS]
        if not sp.issparse(_counts):
            errors.append(f"layers['{LAYER_COUNTS}'] 非稀疏矩阵（got {type(_counts).__name__}）")
        else:
            _data = _counts.data
            if (_data < 0).any():
                errors.append(f"layers['{LAYER_COUNTS}'] 含负值")
            # 近整数校验（允许 float32 浮点误差 < 1e-4）
            _diff = np.abs(_data - np.round(_data))
            if (_diff >= 1e-4).any():
                _max_diff = float(_diff.max())
                errors.append(
                    f"layers['{LAYER_COUNTS}'] 含非整数值（max |x - round(x)| = {_max_diff:.2e}）"
                )

    # 2. doublet_class 列存在 + Categorical + 取值受限
    if DOUBLET_CLASS_COL not in adata.obs.columns:
        errors.append(f"obs['{DOUBLET_CLASS_COL}'] 缺失")
    else:
        _dbl_series = adata.obs[DOUBLET_CLASS_COL]
        if not isinstance(_dbl_series.dtype, pd.CategoricalDtype):
            errors.append(
                f"obs['{DOUBLET_CLASS_COL}'] dtype 应为 Categorical，实际为 {_dbl_series.dtype}"
            )
        else:
            _seen = set(_dbl_series.cat.categories)
            _invalid = _seen - _VALID_DOUBLET_CLASSES
            if _invalid:
                errors.append(
                    f"obs['{DOUBLET_CLASS_COL}'] 含非法取值 {sorted(_invalid)}，"
                    f"合法取值仅限 {sorted(_VALID_DOUBLET_CLASSES)}"
                )
            _present = set(_dbl_series.dropna().unique())
            _unexpected = _present - _VALID_DOUBLET_CLASSES
            if _unexpected:
                errors.append(
                    f"obs['{DOUBLET_CLASS_COL}'] 实际数据含非法取值 {sorted(_unexpected)}"
                )

    # 3. doublet_include 列存在 + bool dtype
    if DOUBLET_INCLUDE_COL not in adata.obs.columns:
        errors.append(f"obs['{DOUBLET_INCLUDE_COL}'] 缺失")
    else:
        if adata.obs[DOUBLET_INCLUDE_COL].dtype != bool:
            errors.append(
                f"obs['{DOUBLET_INCLUDE_COL}'] dtype 应为 bool，实际为 {adata.obs[DOUBLET_INCLUDE_COL].dtype}"
            )

    # 4. expression_contract 在 uns 中
    if EXPRESSION_CONTRACT_KEY not in adata.uns:
        errors.append(f"uns['{EXPRESSION_CONTRACT_KEY}'] 缺失")
    else:
        _ec = adata.uns[EXPRESSION_CONTRACT_KEY]
        if not isinstance(_ec, dict):
            errors.append(
                f"uns['{EXPRESSION_CONTRACT_KEY}'] 应为 dict，实际为 {type(_ec).__name__}"
            )

    # 5. doublet_contract 在 uns 中（若 doublet 检测曾运行）
    if DOUBLET_CLASS_COL in adata.obs.columns:
        _has_non_singlet = (adata.obs[DOUBLET_CLASS_COL] != "singlet").any()
        if _has_non_singlet and DOUBLET_CONTRACT_KEY not in adata.uns:
            warnings.append(
                f"uns['{DOUBLET_CONTRACT_KEY}'] 缺失：存在非 singlet 细胞但未记录 doublet contract"
            )
        elif DOUBLET_CONTRACT_KEY in adata.uns:
            _dc = adata.uns[DOUBLET_CONTRACT_KEY]
            if not isinstance(_dc, dict):
                warnings.append(
                    f"uns['{DOUBLET_CONTRACT_KEY}'] 应为 dict，实际为 {type(_dc).__name__}"
                )

    # 6. adata.var 包含 gene_ids 或 symbol 列（至少一个）
    # guard: mock/synthetic adata may have .var = None
    if hasattr(adata, "var") and adata.var is not None and hasattr(adata.var, "columns"):
        _has_gene_ids = any(
            col.lower() in ("gene_ids", "ensembl_id", "ensembl")
            for col in adata.var.columns
        )
        _has_symbol = any(
            col.lower() in ("symbol", "gene_symbol", "gene_name")
            for col in adata.var.columns
        )
        if not _has_gene_ids and not _has_symbol:
            warnings.append(
                "adata.var 缺少基因标识列（需要 gene_ids/ensembl_id 或 symbol/gene_symbol 中的至少一个）"
            )

    passed = len(errors) == 0
    return {"passed": passed, "errors": errors, "warnings": warnings}
