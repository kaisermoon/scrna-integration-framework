"""Tests for scorer functions — synthetic AnnData, no real data files."""

import anndata
import numpy as np

from scrna_integration.scorers import (
    annotation_concordance,
    clustering_metrics,
    integration_metrics,
    qc_balance,
)


def _make_adata(n_cells=200, n_genes=100, n_pcs=10, seed=42):
    rng = np.random.default_rng(seed)
    X = rng.poisson(2, size=(n_cells, n_genes)).astype(np.float32)  # noqa: N806
    adata = anndata.AnnData(X)
    # Simulate 2-cluster PCA embedding
    pcs = rng.normal(size=(n_cells, n_pcs)).astype(np.float32)
    pcs[: n_cells // 2] += 3.0
    adata.obsm["X_pca"] = pcs
    # QC metrics
    adata.obs["n_genes_by_counts"] = rng.integers(200, 5000, size=n_cells).astype(np.float64)
    return adata


class TestQcBalance:
    def test_returns_retention_ratios(self):
        before = _make_adata(200, 100)
        after = before[:150, :80].copy()
        r = qc_balance(after, before, {})
        assert 0.6 < r["cell_retention"] < 0.8
        assert r["gene_retention"] == 0.8
        assert "mean_genes_after" in r
        assert "mean_genes_before" in r

    def test_empty_before_does_not_divide_by_zero(self):
        before = _make_adata(0)
        after = _make_adata(100)
        r = qc_balance(after, before, {})
        assert r["cell_retention"] > 0
        assert r["gene_retention"] > 0


class TestClusteringMetrics:
    def test_silhouette_with_pca_and_leiden(self):
        adata = _make_adata(200)
        adata.obs["leiden"] = ["0"] * 100 + ["1"] * 100
        r = clustering_metrics(adata, adata, {})
        assert "silhouette" in r
        assert r["silhouette"] > 0.3  # well-separated clusters

    def test_ari_when_known_labels_present(self):
        adata = _make_adata(200)
        adata.obs["leiden"] = ["0"] * 100 + ["1"] * 100
        adata.obs["cell_type"] = ["A"] * 100 + ["B"] * 100
        r = clustering_metrics(adata, adata, {})
        assert "ari" in r
        assert r["ari"] == 1.0

    def test_graceful_on_single_cluster(self):
        adata = _make_adata(200)
        adata.obs["leiden"] = ["0"] * 200
        r = clustering_metrics(adata, adata, {})
        assert "silhouette" not in r
        assert "_note" in r

    def test_no_pca_no_crash(self):
        adata = _make_adata(200)
        del adata.obsm["X_pca"]
        adata.obs["leiden"] = ["0"] * 100 + ["1"] * 100
        r = clustering_metrics(adata, adata, {})
        assert isinstance(r, dict)
        assert "silhouette" not in r

    def test_no_cluster_col_returns_early(self):
        adata = _make_adata(200)
        r = clustering_metrics(adata, adata, {})
        assert "_note" in r


class TestIntegrationMetrics:
    def test_returns_silhouette_metrics(self):
        adata = _make_adata(200)
        adata.obs["batch"] = ["A"] * 100 + ["B"] * 100
        # Interleaved cell types across batches for realistic scenario
        adata.obs["cell_type"] = (["T"] * 50 + ["B"] * 50) * 2
        r = integration_metrics(adata, adata, {})
        assert "silhouette_batch" in r
        assert "silhouette_celltype" in r
        assert "scib_available" in r

    def test_no_batch_column_returns_gracefully(self):
        adata = _make_adata(200)
        r = integration_metrics(adata, adata, {})
        assert "_note" in r

    def test_no_embedding_returns_gracefully(self):
        adata = _make_adata(200)
        del adata.obsm["X_pca"]
        adata.obs["batch"] = ["A"] * 200
        r = integration_metrics(adata, adata, {})
        assert isinstance(r, dict)
        assert "silhouette_batch" not in r

    def test_scib_optional_not_hard_dependency(self):
        """Integration metrics does not crash when scib-metrics is missing."""
        adata = _make_adata(200)
        adata.obs["batch"] = ["A"] * 100 + ["B"] * 100
        r = integration_metrics(adata, adata, {})
        # scib is not installed → scib_available == 0.0
        assert r["scib_available"] == 0.0


class TestAnnotationConcordance:
    def test_kappa_on_two_label_columns(self):
        adata = _make_adata(200)
        adata.obs["label_a"] = ["X"] * 100 + ["Y"] * 100
        adata.obs["label_b"] = ["X"] * 95 + ["Y"] * 5 + ["X"] * 5 + ["Y"] * 95
        r = annotation_concordance(adata, adata, {})
        assert "cohen_kappa" in r
        assert 0.7 < r["cohen_kappa"] < 1.0

    def test_perfect_agreement_kappa_one(self):
        adata = _make_adata(200)
        adata.obs["cell_type_a"] = ["X"] * 100 + ["Y"] * 100
        adata.obs["cell_type_b"] = ["X"] * 100 + ["Y"] * 100
        r = annotation_concordance(adata, adata, {})
        assert r["cohen_kappa"] == 1.0

    def test_missing_columns_returns_gracefully(self):
        adata = _make_adata(200)
        r = annotation_concordance(adata, adata, {})
        assert "_note" in r

    def test_params_specify_columns(self):
        adata = _make_adata(200)
        adata.obs["foo"] = ["A"] * 100 + ["B"] * 100
        adata.obs["bar"] = ["A"] * 95 + ["B"] * 5 + ["A"] * 5 + ["B"] * 95
        r = annotation_concordance(adata, adata, {"label_a": "foo", "label_b": "bar"})
        assert "cohen_kappa" in r
