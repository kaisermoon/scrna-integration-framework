"""scorer 函数测试——合成 AnnData，无需真实数据文件。

签名为 ADR-0009 后的直接调用形式：单 AnnData + 可选 key 参数。
"""

import anndata
import numpy as np

from scrna_integration.scorers import (
    annotation_concordance,
    clustering_metrics,
    integration_metrics,
    qc_balance,
)


def _make_adata(n_cells=200, n_genes=100, n_pcs=10, seed=42):
    """构建测试用合成 AnnData。"""
    rng = np.random.default_rng(seed)
    X = rng.poisson(2, size=(n_cells, n_genes)).astype(np.float32)  # noqa: N806
    adata = anndata.AnnData(X)
    # 模拟 2 簇 PCA 嵌入
    pcs = rng.normal(size=(n_cells, n_pcs)).astype(np.float32)
    pcs[: n_cells // 2] += 3.0
    adata.obsm["X_pca"] = pcs
    # QC 指标
    adata.obs["n_genes_by_counts"] = rng.integers(200, 5000, size=n_cells).astype(np.float64)
    return adata


class TestQcBalance:
    """QC 平衡指标：过滤前后细胞/基因保留比例。"""

    def test_returns_retention_ratios(self):
        before = _make_adata(200, 100)
        after = before[:150, :80].copy()
        r = qc_balance(after, before)
        assert 0.6 < r["cell_retention"] < 0.8
        assert r["gene_retention"] == 0.8
        assert "mean_genes_after" in r
        assert "mean_genes_before" in r

    def test_empty_before_does_not_divide_by_zero(self):
        before = _make_adata(0)
        after = _make_adata(100)
        r = qc_balance(after, before)
        assert r["cell_retention"] > 0
        assert r["gene_retention"] > 0


class TestClusteringMetrics:
    """聚类质量指标：Silhouette + ARI。"""

    def test_silhouette_with_pca_and_leiden(self):
        adata = _make_adata(200)
        adata.obs["leiden"] = ["0"] * 100 + ["1"] * 100
        r = clustering_metrics(adata)
        assert "silhouette" in r
        assert r["silhouette"] > 0.3  # 分离良好的簇

    def test_ari_when_known_labels_present(self):
        adata = _make_adata(200)
        adata.obs["leiden"] = ["0"] * 100 + ["1"] * 100
        adata.obs["cell_type"] = ["A"] * 100 + ["B"] * 100
        r = clustering_metrics(adata)
        assert "ari" in r
        assert r["ari"] == 1.0

    def test_graceful_on_single_cluster(self):
        adata = _make_adata(200)
        adata.obs["leiden"] = ["0"] * 200
        r = clustering_metrics(adata)
        assert "silhouette" not in r
        assert "_note" in r

    def test_no_pca_no_crash(self):
        adata = _make_adata(200)
        del adata.obsm["X_pca"]
        adata.obs["leiden"] = ["0"] * 100 + ["1"] * 100
        r = clustering_metrics(adata)
        assert isinstance(r, dict)
        assert "silhouette" not in r

    def test_no_cluster_col_returns_early(self):
        adata = _make_adata(200)
        r = clustering_metrics(adata)
        assert "_note" in r

    def test_explicit_cluster_key(self):
        """指定 cluster_key 参数可显式指定聚类列。"""
        adata = _make_adata(200)
        adata.obs["my_clusters"] = ["A"] * 100 + ["B"] * 100
        r = clustering_metrics(adata, cluster_key="my_clusters")
        assert "silhouette" in r

    def test_explicit_label_key(self):
        """指定 label_key 参数可显式指定已知标签列。"""
        adata = _make_adata(200)
        adata.obs["leiden"] = ["0"] * 100 + ["1"] * 100
        adata.obs["ground_truth"] = ["X"] * 100 + ["Y"] * 100
        r = clustering_metrics(adata, label_key="ground_truth")
        assert "ari" in r
        assert r["ari"] == 1.0


class TestIntegrationMetrics:
    """批次整合指标：batch 混合 + celltype 保留。"""

    def test_returns_silhouette_metrics(self):
        adata = _make_adata(200)
        adata.obs["batch"] = ["A"] * 100 + ["B"] * 100
        # 批次间交错排列细胞类型以模拟真实场景
        adata.obs["cell_type"] = (["T"] * 50 + ["B"] * 50) * 2
        r = integration_metrics(adata)
        assert "silhouette_batch" in r
        assert "silhouette_celltype" in r
        assert "scib_available" in r

    def test_no_batch_column_returns_gracefully(self):
        adata = _make_adata(200)
        r = integration_metrics(adata, batch_key="nonexistent_column")
        assert "_note" in r

    def test_no_embedding_returns_gracefully(self):
        adata = _make_adata(200)
        del adata.obsm["X_pca"]
        adata.obs["batch"] = ["A"] * 200
        r = integration_metrics(adata)
        assert isinstance(r, dict)
        assert "silhouette_batch" not in r

    def test_scib_optional_not_hard_dependency(self):
        """scib-metrics 缺失时 integration_metrics 不崩溃。"""
        adata = _make_adata(200)
        adata.obs["batch"] = ["A"] * 100 + ["B"] * 100
        r = integration_metrics(adata)
        # scib 未安装 → scib_available == 0.0
        assert r["scib_available"] == 0.0

    def test_explicit_label_key(self):
        """显式指定 label_key 绕过自动检测。"""
        adata = _make_adata(200)
        adata.obs["batch"] = ["A"] * 100 + ["B"] * 100
        adata.obs["custom_ct"] = (["T"] * 50 + ["B"] * 50) * 2
        r = integration_metrics(adata, label_key="custom_ct")
        assert "silhouette_celltype" in r


class TestAnnotationConcordance:
    """标注一致性：Cohen's kappa。"""

    def test_kappa_on_two_label_columns(self):
        adata = _make_adata(200)
        adata.obs["label_a"] = ["X"] * 100 + ["Y"] * 100
        adata.obs["label_b"] = ["X"] * 95 + ["Y"] * 5 + ["X"] * 5 + ["Y"] * 95
        r = annotation_concordance(adata)
        assert "cohen_kappa" in r
        assert 0.7 < r["cohen_kappa"] < 1.0

    def test_perfect_agreement_kappa_one(self):
        adata = _make_adata(200)
        adata.obs["cell_type_a"] = ["X"] * 100 + ["Y"] * 100
        adata.obs["cell_type_b"] = ["X"] * 100 + ["Y"] * 100
        r = annotation_concordance(adata)
        assert r["cohen_kappa"] == 1.0

    def test_missing_columns_returns_gracefully(self):
        adata = _make_adata(200)
        r = annotation_concordance(adata)
        assert "_note" in r

    def test_explicit_label_columns(self):
        adata = _make_adata(200)
        adata.obs["foo"] = ["A"] * 100 + ["B"] * 100
        adata.obs["bar"] = ["A"] * 95 + ["B"] * 5 + ["A"] * 5 + ["B"] * 95
        r = annotation_concordance(adata, label_a="foo", label_b="bar")
        assert "cohen_kappa" in r
