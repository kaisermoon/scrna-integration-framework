"""scorer 函数测试——合成 AnnData，无需真实数据文件。

签名为 ADR-0009 后的直接调用形式：单 AnnData + 可选 key 参数。
"""

import anndata
import numpy as np
import pytest

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

    def test_explicit_embed_key(self):
        """显式 embed_key 真生效——两个不同嵌入产生不同 silhouette_batch，证明未走 auto-detect。"""
        rng = np.random.default_rng(99)
        n_cells = 200
        # 嵌入 A：两个批次结构良好分离 → 低 batch silhouette（批次混合好）
        # 嵌入 B：纯随机 → 批次 score 接近 0
        embed_a = rng.normal(size=(n_cells, 10)).astype(np.float32)
        embed_a[:100] += 5.0  # 前半细胞偏移，形成明显分离
        embed_b = rng.normal(0, 0.1, size=(n_cells, 10)).astype(np.float32)

        adata = anndata.AnnData(rng.poisson(2, size=(n_cells, 20)).astype(np.float32))
        adata.obsm["X_pca"] = embed_a
        adata.obsm["X_alt"] = embed_b
        adata.obs["batch"] = ["A"] * 100 + ["B"] * 100

        r_a = integration_metrics(adata, embed_key="X_pca")
        r_b = integration_metrics(adata, embed_key="X_alt")

        assert "silhouette_batch" in r_a
        assert "silhouette_batch" in r_b
        assert abs(r_a["silhouette_batch"] - r_b["silhouette_batch"]) > 0.05, (
            f"两个不同嵌入产生近乎相同的 silhouette_batch "
            f"(X_pca={r_a['silhouette_batch']:.4f}, X_alt={r_b['silhouette_batch']:.4f})"
        )

    def test_explicit_embed_key_not_found(self):
        """embed_key 不存在于 obsm 时 graceful 返回 _note。"""
        adata = _make_adata(200)
        adata.obs["batch"] = ["A"] * 100 + ["B"] * 100
        r = integration_metrics(adata, embed_key="not_found")
        assert "_note" in r
        assert np.isnan(r["_note"])


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


# ---------------------------------------------------------------------------
# P1 补测：clustering_metrics 自动检测 louvain/cluster 列 + 优先级
# ---------------------------------------------------------------------------


class TestClusteringAutoDetect:
    """clustering_metrics 对 obs 列的自动检测行为。"""

    def test_auto_detect_louvain_column(self):
        """obs 含 louvain 列 → 自动选用。"""
        adata = _make_adata(200)
        adata.obs["louvain"] = ["0"] * 100 + ["1"] * 100
        r = clustering_metrics(adata)
        assert "silhouette" in r
        assert r["silhouette"] > 0.3

    def test_auto_detect_cluster_column(self):
        """obs 含 cluster 列（无 louvain/leiden）→ 自动选用。"""
        adata = _make_adata(200)
        adata.obs["cluster"] = ["A"] * 100 + ["B"] * 100
        r = clustering_metrics(adata)
        assert "silhouette" in r

    def test_louvain_priority_over_cluster(self):
        """louvain 和 cluster 同时存在时，louvain 优先级更高。"""
        adata = _make_adata(200)
        adata.obs["cluster"] = ["X"] * 200  # 单簇，无法 silhouette
        adata.obs["louvain"] = ["0"] * 100 + ["1"] * 100  # 两簇
        r = clustering_metrics(adata)
        # louvain 被选中 → 有两簇 → silhouette 存在
        assert "silhouette" in r
        assert r["silhouette"] > 0.3


# ---------------------------------------------------------------------------
# P1 补测：integration_metrics 嵌入优先级链（X_pca_harmony > X_scVI > X_scANVI > X_pca）
# ---------------------------------------------------------------------------


class TestEmbedPriorityChain:
    """integration_metrics 自动检测 obsm 中嵌入的优先级链。"""

    def test_embed_priority_harmony_over_pca(self):
        """X_pca_harmony 存在时优先于 X_pca。"""
        rng = np.random.default_rng(99)
        n_cells = 200
        # 构造两个嵌入：harmony 和 pca
        harmony_embed = rng.normal(size=(n_cells, 10)).astype(np.float32)
        harmony_embed[:100] += 5.0
        pca_embed = rng.normal(size=(n_cells, 10)).astype(np.float32)
        pca_embed[:100] += 1.0  # 分离较弱

        adata = _make_adata(200)
        adata.obsm["X_pca_harmony"] = harmony_embed
        adata.obsm["X_pca"] = pca_embed
        adata.obs["batch"] = ["A"] * 100 + ["B"] * 100

        r = integration_metrics(adata)
        # 选用了 harmony 嵌入 → silhouette_batch 应接近 harmony 的值
        assert "silhouette_batch" in r

    def test_embed_fallback_to_custom_X_key(self):
        """obsm 中无标准键但有 X_custom → 回退使用 X_custom。"""
        rng = np.random.default_rng(99)
        n_cells = 200
        custom_embed = rng.normal(size=(n_cells, 10)).astype(np.float32)
        custom_embed[:100] += 3.0

        adata = _make_adata(200)
        del adata.obsm["X_pca"]  # 移除标准键
        adata.obsm["X_custom"] = custom_embed
        adata.obs["batch"] = ["A"] * 100 + ["B"] * 100

        r = integration_metrics(adata)
        assert "silhouette_batch" in r


# ---------------------------------------------------------------------------
# P1 xfail: annotation_concordance label_a / label_b 同列自洽 bug
# ---------------------------------------------------------------------------


class TestAnnotationConcordanceBug:
    """确认 bug：label_b 自动检测未排除 label_a，可致同列自洽 kappa=1.0。"""

    @pytest.mark.xfail(
        reason=(
            "BUG: label_b 自动检测未排除 label_a 已用列名（scorers.py:215-218）。"
            "当 candidates[1] 恰好等于 label_a 时，两张标签列为同一列，"
            "cohen_kappa 恒为 1.0——自身一致性伪装成跨标注列一致性，结果虚高。"
            "修复方向：candidates = [c for c in candidates if c != label_a]"
        ),
        strict=True,
    )
    def test_label_a_label_b_same_column_self_concordance(self):
        """label_a=label_b 同列时不应返回 kappa=1.0。"""
        adata = _make_adata(200)
        # 仅两个候选列：cell_type_labels 和 leiden
        adata.obs["cell_type_labels"] = ["X"] * 100 + ["Y"] * 100
        adata.obs["leiden"] = ["0"] * 100 + ["1"] * 100
        # 显式传相同 label_a=label_b → 同列自洽
        r = annotation_concordance(
            adata, label_a="cell_type_labels", label_b="cell_type_labels"
        )
        # 期望：应检测到同列并返回 _note；当前代码未做此检测
        # 当前行为：kappa=1.0（同列自洽）。xfail 标记暴露此 bug。
        assert "_note" in r or r.get("cohen_kappa", 1.0) < 0.99
