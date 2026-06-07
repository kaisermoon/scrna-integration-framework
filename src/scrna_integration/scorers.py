"""Scorer functions for the sweep harness.

Each scorer has the signature ``(adata_after, adata_before, params) -> dict[str, float]``.
They are plain functions, not plugin classes (ADR-0003).
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


def qc_balance(
    adata_after: AnnData, adata_before: AnnData, params: dict
) -> dict[str, float]:
    """Cell and gene retention ratios after QC filtering."""
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
    adata_after: AnnData, adata_before: AnnData, params: dict
) -> dict[str, float]:
    """Silhouette score (on PCA) + ARI if known cell-type labels exist.

    Returns gracefully (empty or partial dict) when cluster labels are
    absent, single-valued, or no PCA embedding is available.
    """
    result: dict[str, float] = {}

    # Locate cluster labels (scanpy default naming)
    cluster_col = _first_match(adata_after.obs.columns, ["leiden", "louvain", "cluster"])
    if cluster_col is None:
        result["_note"] = float("nan")
        return result

    labels = adata_after.obs[cluster_col].astype(str)
    mask = adata_after.obs[cluster_col].notna()
    valid_labels = labels[mask]

    if len(np.unique(valid_labels)) < 2:
        result["_note"] = float("nan")
        return result

    # Silhouette on PCA embedding
    if "X_pca" in adata_after.obsm:
        x = adata_after[mask].obsm["X_pca"]  # noqa: N806
        try:
            result["silhouette"] = float(silhouette_score(x, valid_labels))
        except Exception:
            pass

    # ARI against known labels if available
    true_col = _first_match(adata_after.obs.columns, ["cell_type", "label", "annotation"])
    if true_col is not None:
        true_labels = adata_after.obs[true_col][mask].astype(str)
        try:
            result["ari"] = float(adjusted_rand_score(true_labels, valid_labels))
        except Exception:
            pass

    return result


def integration_metrics(
    adata_after: AnnData, adata_before: AnnData, params: dict
) -> dict[str, float]:
    """Batch-mixing metrics: silhouette by batch, silhouette by cell-type.

    ``scib-metrics`` is optional — when absent only the sklearn subset is
    computed and ``scib_available`` is set to 0.0.
    """
    result: dict[str, float] = {}

    batch_col = _first_match(adata_after.obs.columns, ["batch", "sample_id", "source_dataset"])
    if batch_col is None:
        return {"_note": float("nan")}

    # Locate any embedding (prefer integration output, fallback PCA)
    embed_key = _first_match(
        list(adata_after.obsm.keys()),
        ["X_pca_harmony", "X_scVI", "X_scANVI", "X_pca"],
    )
    if embed_key is None:
        # Try any obsm key starting with X_
        for k in adata_after.obsm.keys():
            if k.startswith("X_"):
                embed_key = k
                break
    if embed_key is None:
        return {"_note": float("nan")}

    x_embed = adata_after.obsm[embed_key]  # noqa: N806

    # Silhouette by batch (lower = better mixing)
    batch_labels = adata_after.obs[batch_col].astype(str)
    bm = adata_after.obs[batch_col].notna()
    if bm.sum() >= 2 and len(np.unique(batch_labels[bm])) >= 2:
        try:
            result["silhouette_batch"] = float(silhouette_score(x_embed[bm], batch_labels[bm]))
        except Exception:
            pass

    # Silhouette by cell-type (higher = better biology preservation)
    ct_col = _first_match(
        adata_after.obs.columns,
        ["cell_type", "cell_type_final_v1", "label", "cell_type_original"],
    )
    if ct_col is not None:
        ct_labels = adata_after.obs[ct_col].astype(str)
        cm = adata_after.obs[ct_col].notna()
        if cm.sum() >= 2 and len(np.unique(ct_labels[cm])) >= 2:
            try:
                result["silhouette_celltype"] = float(
                    silhouette_score(x_embed[cm], ct_labels[cm])
                )
            except Exception:
                pass

    # scib-metrics availability marker
    try:
        import scib_metrics  # noqa: F401

        result["scib_available"] = 1.0
    except ImportError:
        result["scib_available"] = 0.0

    return result


def annotation_concordance(
    adata_after: AnnData, adata_before: AnnData, params: dict
) -> dict[str, float]:
    """Cohen's kappa between two label columns in ``adata_after.obs``.

    Columns are taken from *params* (``label_a`` / ``label_b``) or
    auto-detected from common naming patterns.
    """
    col_a = params.get("label_a")
    col_b = params.get("label_b")

    if col_a is None or col_b is None:
        # Auto-detect two label-like columns
        candidates = [
            c
            for c in adata_after.obs.columns
            if any(p in c.lower() for p in ("cell_type", "label", "annotation", "leiden"))
        ]
        if len(candidates) >= 2 and col_a is None:
            col_a = candidates[0]
        if len(candidates) >= 2 and col_b is None:
            col_b = candidates[1]

    if col_a is None or col_b is None or col_a not in adata_after.obs.columns or col_b not in adata_after.obs.columns:
        return {"_note": float("nan")}

    a = adata_after.obs[col_a].astype(str)
    b = adata_after.obs[col_b].astype(str)
    valid = adata_after.obs[col_a].notna() & adata_after.obs[col_b].notna()

    return {"cohen_kappa": float(cohen_kappa_score(a[valid], b[valid]))}


def _first_match(columns: list[str], candidates: list[str]) -> str | None:
    """Return the first column in *columns* whose name contains any *candidate*."""
    for candidate in candidates:
        for col in columns:
            if candidate in col:
                return col
    return None
