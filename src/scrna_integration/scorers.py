"""可直接调用的指标函数（无回调抽象，面向非计算机专业 PI/学生）。

每个函数接收一个 AnnData + 可选的 key 参数，返回 ``dict[str, float]``。
在 notebook 的显式 for 循环中直接调用：
   m = integration_metrics(adata_copy, batch_key="batch")
   results.append({"use_rep": rep, **m})

与 ADR-0003 一致：纯函数，无 plugin 类文件，无注册中心。
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from sklearn.metrics import (
    adjusted_rand_score,
    cohen_kappa_score,
    silhouette_score,
)

if TYPE_CHECKING:
    from anndata import AnnData


def qc_balance(adata_after: AnnData, adata_before: AnnData) -> dict[str, float]:
    """QC 过滤前后细胞数与基因数保留比例。

    Args:
        adata_after: 过滤后的 AnnData
        adata_before: 过滤前的 AnnData

    Returns:
        cell_retention / gene_retention / mean_genes_after / mean_genes_before
    """
    result: dict[str, float] = {
        "cell_retention": adata_after.n_obs / max(adata_before.n_obs, 1),
        "gene_retention": adata_after.n_vars / max(adata_before.n_vars, 1),
    }
    if "n_genes_by_counts" in adata_after.obs.columns:
        result["mean_genes_after"] = float(adata_after.obs["n_genes_by_counts"].mean())
    if "n_genes_by_counts" in adata_before.obs.columns:
        result["mean_genes_before"] = float(adata_before.obs["n_genes_by_counts"].mean())
    return result


def clustering_metrics(
    adata: AnnData,
    cluster_key: str | None = None,
    label_key: str | None = None,
) -> dict[str, float]:
    """聚类质量指标：轮廓系数（Silhouette）+ 调整兰德指数（ARI）。

    自动检测 obs 中的 leiden/louvain/cluster 列作为聚类标签，
    以及 cell_type/label/annotation 列作为已知标签计算 ARI。

    Args:
        adata: 含 obsm['X_pca'] 和 obs 聚类列的 AnnData
        cluster_key: 指定聚类标签列名（None 则自动检测）
        label_key: 指定已知标签列名（None 则自动检测）

    Returns:
        silhouette（PCA 空间轮廓系数）/ ari（如已知标签存在）
        无可聚类列时返回 {"_note": nan}
    """
    result: dict[str, float] = {}

    # 定位聚类标签列（scanpy 默认命名）
    if cluster_key is None:
        cluster_key = _first_match(adata.obs.columns, ["leiden", "louvain", "cluster"])
    if cluster_key is None:
        result["_note"] = float("nan")
        return result

    labels = adata.obs[cluster_key].astype(str)
    mask = adata.obs[cluster_key].notna()
    valid_labels = labels[mask]

    if len(np.unique(valid_labels)) < 2:
        result["_note"] = float("nan")
        return result

    # PCA 空间的轮廓系数（Silhouette）
    if "X_pca" in adata.obsm:
        x = adata[mask].obsm["X_pca"]  # noqa: N806
        try:
            result["silhouette"] = float(silhouette_score(x, valid_labels))
        except Exception:
            pass

    # ARI：与已知标签比较
    if label_key is None:
        label_key = _first_match(adata.obs.columns, ["cell_type", "label", "annotation"])
    if label_key is not None:
        true_labels = adata.obs[label_key][mask].astype(str)
        try:
            result["ari"] = float(adjusted_rand_score(true_labels, valid_labels))
        except Exception:
            pass

    return result


def integration_metrics(
    adata: AnnData,
    batch_key: str = "batch",
    label_key: str | None = None,
    embed_key: str | None = None,
) -> dict[str, float]:
    """批次整合指标：按批次和细胞类型计算轮廓系数，评估批次混合与生物学信号保留。

    当 *embed_key* 显式给出时直接使用该嵌入；否则自动检测 obsm 中的嵌入
    （优先 X_pca_harmony/X_scVI/X_scANVI，回退到任意 X_* 键）。

    ``scib-metrics`` 可选——未安装时仅计算 sklearn 子集并将 ``scib_available`` 设为 0.0。

    Args:
        adata: 含 obsm 嵌入和 obs 批次/标签列的 AnnData
        batch_key: obs 批次列名（默认 "batch"）
        label_key: obs 细胞类型列名（None 则自动检测 cell_type/cell_type_final_v1）
        embed_key: 嵌入 obsm 键名（None 则按优先级自动检测）

    Returns:
        silhouette_batch（越低=批次混合越好）/
        silhouette_celltype（越高=生物学信号保留越好）/
        scib_available（1.0=scib-metrics 已安装，0.0=未安装）
    """
    result: dict[str, float] = {}

    batch_col = _first_match(adata.obs.columns, [batch_key, "sample_id", "source_dataset"])
    if batch_col is None:
        return {"_note": float("nan")}

    # 定位嵌入：优先使用显式给定的 embed_key，否则自动检测
    if embed_key is not None:
        if embed_key not in adata.obsm:
            return {"_note": float("nan")}
    else:
        embed_key = _first_match(
            list(adata.obsm.keys()),
            ["X_pca_harmony", "X_scVI", "X_scANVI", "X_pca"],
        )
        if embed_key is None:
            # 回退：任意 X_ 开头的 obsm 键
            for k in adata.obsm.keys():
                if k.startswith("X_"):
                    embed_key = k
                    break
        if embed_key is None:
            return {"_note": float("nan")}

    x_embed = adata.obsm[embed_key]  # noqa: N806

    # 按批次轮廓系数（越低=批次混合越好）
    batch_labels = adata.obs[batch_col].astype(str)
    bm = adata.obs[batch_col].notna()
    if bm.sum() >= 2 and len(np.unique(batch_labels[bm])) >= 2:
        try:
            result["silhouette_batch"] = float(silhouette_score(x_embed[bm], batch_labels[bm]))
        except Exception:
            pass

    # 按细胞类型轮廓系数（越高=生物学信号保留越好）
    if label_key is None:
        label_key = _first_match(
            adata.obs.columns,
            ["cell_type", "cell_type_final_v1", "label", "cell_type_original"],
        )
    if label_key is not None:
        ct_labels = adata.obs[label_key].astype(str)
        cm = adata.obs[label_key].notna()
        if cm.sum() >= 2 and len(np.unique(ct_labels[cm])) >= 2:
            try:
                result["silhouette_celltype"] = float(
                    silhouette_score(x_embed[cm], ct_labels[cm])
                )
            except Exception:
                pass

    # scib-metrics 可用性标记
    try:
        import scib_metrics  # noqa: F401

        result["scib_available"] = 1.0
    except ImportError:
        result["scib_available"] = 0.0

    return result


def annotation_concordance(
    adata: AnnData,
    label_a: str | None = None,
    label_b: str | None = None,
) -> dict[str, float]:
    """两个标注列之间的 Cohen's kappa 一致性系数。

    自动检测 obs 中包含 cell_type/label/annotation/leiden 的列作为两个标注列。

    Args:
        adata: 含两个标注列的 AnnData
        label_a: 第一个标注列名（None 则自动检测）
        label_b: 第二个标注列名（None 则自动检测）

    Returns:
        cohen_kappa（-1 到 1，1=完美一致，0=偶然一致）
    """
    if label_a is None or label_b is None:
        # 自动检测两个类似标注的列
        candidates = [
            c
            for c in adata.obs.columns
            if any(p in c.lower() for p in ("cell_type", "label", "annotation", "leiden"))
        ]
        if len(candidates) >= 2 and label_a is None:
            label_a = candidates[0]
        if len(candidates) >= 2 and label_b is None:
            label_b = candidates[1]

    if label_a is None or label_b is None or label_a not in adata.obs.columns or label_b not in adata.obs.columns:
        return {"_note": float("nan")}

    a = adata.obs[label_a].astype(str)
    b = adata.obs[label_b].astype(str)
    valid = adata.obs[label_a].notna() & adata.obs[label_b].notna()

    return {"cohen_kappa": float(cohen_kappa_score(a[valid], b[valid]))}


def _first_match(columns: list[str], candidates: list[str]) -> str | None:
    """返回 *columns* 中第一个名称包含 *candidates* 中任一项的列。"""
    for candidate in candidates:
        for col in columns:
            if candidate in col:
                return col
    return None
