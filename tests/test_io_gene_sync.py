"""_is_ensembl_id / _infer_species_from_ensembl_prefix /
_infer_species_from_batch_ensembl_ids / _sync_ensembl_to_symbol 的单元测试。

mygene 查询全部 mock，绝不允许测试中发起网络请求。
"""

from __future__ import annotations

import warnings
from unittest import mock

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp

from scrna_integration.io import (
    _infer_species_from_batch_ensembl_ids,
    _infer_species_from_ensembl_prefix,
    _is_ensembl_id,
    _sync_ensembl_to_symbol,
    sync_gene_ids,
)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _make_adata(var_index, var_columns=None, uns=None):
    """构造最小 AnnData，含有简化的 X 矩阵和 obs。"""
    n_genes = len(var_index)
    var = pd.DataFrame(index=var_index)
    if var_columns:
        for col_name, col_values in var_columns.items():
            var[col_name] = col_values

    # 构造最小的整数 counts 矩阵（poisson 采样一小簇细胞）
    rng = np.random.default_rng(42)
    n_cells = 10
    X = rng.poisson(5, size=(n_cells, n_genes)).astype(np.float32)  # noqa: N806
    X = sp.csr_matrix(X)  # noqa: N806

    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    if uns:
        adata.uns.update(uns)
    return adata


# ---------------------------------------------------------------------------
# _is_ensembl_id
# ---------------------------------------------------------------------------


class TestIsEnsemblId:
    def test_ensg_human(self):
        assert _is_ensembl_id("ENSG00000141510") is True

    def test_ensmusg_mouse(self):
        assert _is_ensembl_id("ENSMUSG00000000001") is True

    def test_ensrnog_rat(self):
        assert _is_ensembl_id("ENSRNOG00000000001") is True

    def test_ensdarg_zebrafish(self):
        assert _is_ensembl_id("ENSDARG00000000001") is True

    def test_lowercase_ensg(self):
        assert _is_ensembl_id("ensg00000141510") is True

    def test_non_ensembl_symbol(self):
        assert _is_ensembl_id("TP53") is False

    def test_non_ensembl_gene(self):
        assert _is_ensembl_id("GENE_123") is False

    def test_empty_string(self):
        assert _is_ensembl_id("") is False

    def test_numeric_gene(self):
        assert _is_ensembl_id("12345") is False


# ---------------------------------------------------------------------------
# _infer_species_from_ensembl_prefix
# ---------------------------------------------------------------------------


class TestInferSpeciesFromEnsemblPrefix:
    def test_ensg_human(self):
        assert _infer_species_from_ensembl_prefix("ENSG00000141510") == "human"

    def test_ensmusg_mouse(self):
        assert _infer_species_from_ensembl_prefix("ENSMUSG00000000001") == "mouse"

    def test_ensrnog_rat(self):
        assert _infer_species_from_ensembl_prefix("ENSRNOG00000000001") == "rat"

    def test_unknown_prefix_returns_all(self):
        assert _infer_species_from_ensembl_prefix("ENSXXX00000000001") == "all"

    def test_non_ensembl_returns_all(self):
        assert _infer_species_from_ensembl_prefix("TP53") == "all"

    def test_lowercase_ensmusg(self):
        assert _infer_species_from_ensembl_prefix("ensmusg00000000001") == "mouse"


# ---------------------------------------------------------------------------
# _infer_species_from_batch_ensembl_ids
# ---------------------------------------------------------------------------


class TestInferSpeciesFromBatchEnsemblIds:
    def test_majority_human(self):
        ids = ["ENSG00000000001"] * 15 + ["ENSMUSG00000000001"] * 3
        result = _infer_species_from_batch_ensembl_ids(ids)
        assert result == "human"

    def test_majority_mouse(self):
        ids = ["ENSMUSG00000000001"] * 12 + ["ENSG00000000001"] * 2
        result = _infer_species_from_batch_ensembl_ids(ids)
        assert result == "mouse"

    def test_all_non_ensembl_falls_back_to_adata_uns_species(self):
        adata = _make_adata(["GeneA", "GeneB"], uns={"species": "mouse"})
        result = _infer_species_from_batch_ensembl_ids(
            ["GeneA", "GeneB"], adata=adata
        )
        assert result == "mouse"

    def test_all_non_ensembl_homo_sapiens_normalization(self):
        adata = _make_adata(["GeneA"], uns={"species": "homo_sapiens"})
        result = _infer_species_from_batch_ensembl_ids(["GeneA"], adata=adata)
        assert result == "human"

    def test_all_non_ensembl_mus_musculus_normalization(self):
        adata = _make_adata(["GeneA"], uns={"species": "Mus_musculus"})
        result = _infer_species_from_batch_ensembl_ids(["GeneA"], adata=adata)
        assert result == "mouse"

    def test_all_non_ensembl_rattus_norvegicus_normalization(self):
        adata = _make_adata(["GeneA"], uns={"species": "Rattus_norvegicus"})
        result = _infer_species_from_batch_ensembl_ids(["GeneA"], adata=adata)
        assert result == "rat"

    def test_all_non_ensembl_no_adata_falls_back_human(self):
        result = _infer_species_from_batch_ensembl_ids(["GeneA", "GeneB"])
        assert result == "human"

    def test_empty_list_falls_back_human(self):
        result = _infer_species_from_batch_ensembl_ids([])
        assert result == "human"


# ---------------------------------------------------------------------------
# _sync_ensembl_to_symbol — feature_name 路径（mock mygene，无网络请求）
# ---------------------------------------------------------------------------


class TestSyncEnsemblToSymbolFeatureName:
    """feature_name 列存在且有有效 symbol 的测试。"""

    def test_feature_name_converts_and_preserves_ensembl_id(self):
        """feature_name 中的 symbol 替换 index，原始 Ensembl ID 进 ensembl_id 列。

        使用 mock mygene 防止真实 API 将测试用的 Ensembl ID 映射为实际 symbol。
        """
        mock_client = mock.MagicMock()
        # 二次补全查询只返回 notfound，确保仍为 ENSG 的基因保留不变
        mock_client.querymany.return_value = [
            {"query": "ENSG00000000003", "notfound": True},
        ]

        with mock.patch.dict("sys.modules", {
            "mygene": mock.MagicMock(MyGeneInfo=mock.MagicMock(return_value=mock_client)),
        }):
            adata = _make_adata(
                ["ENSG00000000001", "ENSG00000000002", "ENSG00000000003"],
                var_columns={
                    "feature_name": ["TP53", "BRCA1", "ENSG00000000003"],
                },
            )
            _sync_ensembl_to_symbol(adata)

        assert adata.var.index[0] == "TP53"
        assert adata.var.index[1] == "BRCA1"
        # 无 symbol 的基因保留 Ensembl ID（mygene 二次补全未命中）
        assert adata.var.index[2] == "ENSG00000000003"
        assert adata.var["ensembl_id"].iloc[0] == "ENSG00000000001"
        assert adata.var["ensembl_id"].iloc[1] == "ENSG00000000002"
        assert adata.var["ensembl_id"].iloc[2] == "ENSG00000000003"

    def test_feature_name_with_duplicate_symbols_dedup(self):
        """同一 symbol 出现在多条 feature_name → 去重后 var.index 唯一。"""
        adata = _make_adata(
            ["ENSG00000000001", "ENSG00000000002", "ENSG00000000003"],
            var_columns={
                "feature_name": ["DUPE", "DUPE", "BRCA1"],
            },
        )
        _sync_ensembl_to_symbol(adata)

        assert adata.var.index.is_unique
        # 前两个 DUPE 被去重
        assert adata.var.index[0] == "DUPE"
        assert adata.var.index[1].startswith("DUPE-")  # 去重后缀
        assert adata.var.index[2] == "BRCA1"

    def test_feature_name_path_returns_early_and_skips_mygene_when_all_mapped(self):
        """feature_name 全部是合法 symbol → 不走 mygene 全量兜底。"""
        mock_mg_module = mock.MagicMock()
        mock_mg_module.MyGeneInfo = mock.MagicMock()

        adata = _make_adata(
            ["ENSG00000000001", "ENSG00000000002"],
            var_columns={
                "feature_name": ["TP53", "BRCA1"],
            },
        )
        with mock.patch.dict("sys.modules", {"mygene": mock_mg_module}):
            _sync_ensembl_to_symbol(adata)
            # mygene 全量兜底路径不应该被调用（feature_name 已全覆盖）
            # 且二次补全路径也不应触发（无残留 ENSG）
            mock_mg_module.MyGeneInfo.assert_not_called()

        assert adata.var.index[0] == "TP53"
        assert adata.var.index[1] == "BRCA1"


# ---------------------------------------------------------------------------
# _sync_ensembl_to_symbol — mygene 二次补全（mock）
# ---------------------------------------------------------------------------


class TestSyncEnsemblToSymbolMygeneCompletion:
    """feature_name 中有少量 ENSG 残留的场景 —— mygene 二次补全。"""

    def test_mygene_completion_maps_remaining_ensg(self):
        """feature_name 含 2 个符号 + 1 个 ENSG → mygene 补全那 1 个。"""
        mock_client = mock.MagicMock()
        mock_client.querymany.return_value = [
            {"query": "ENSG00000999999", "symbol": "FAKE"},
        ]

        with mock.patch.dict("sys.modules", {
            "mygene": mock.MagicMock(MyGeneInfo=mock.MagicMock(return_value=mock_client)),
        }):
            adata = _make_adata(
                ["ENSG00000000001", "ENSG00000000002", "ENSG00000999999"],
                var_columns={
                    "feature_name": ["TP53", "BRCA1", "ENSG00000999999"],
                },
            )
            _sync_ensembl_to_symbol(adata)

        assert adata.var.index[0] == "TP53"
        assert adata.var.index[1] == "BRCA1"
        # ENSG00000999999 被 mygene 补全为 FAKE
        assert adata.var.index[2] == "FAKE"

    def test_mygene_completion_handles_all_notfound(self):
        """mygene 二次补全全部未命中 → 保留 ENSG。"""
        mock_client = mock.MagicMock()
        mock_client.querymany.return_value = [
            {"query": "ENSG00000999999", "notfound": True},
        ]

        with mock.patch.dict("sys.modules", {
            "mygene": mock.MagicMock(MyGeneInfo=mock.MagicMock(return_value=mock_client)),
        }):
            adata = _make_adata(
                ["ENSG00000000001", "ENSG00000999999"],
                var_columns={
                    "feature_name": ["TP53", "ENSG00000999999"],
                },
            )
            _sync_ensembl_to_symbol(adata)

        assert adata.var.index[0] == "TP53"
        # ENSG00000999999 保留
        assert adata.var.index[1] == "ENSG00000999999"

    def test_mygene_completion_strips_make_unique_suffix(self):
        """查询前去掉 var_names_make_unique 追加的 -N 后缀。

        模拟场景：ENSG00000999999 同时出现两次，feature_name 去重后变成
        ENSG00000999999 和 ENSG00000999999-1，查询前应分别还原为
        ENSG00000999999 再查 mygene。
        """
        mock_client = mock.MagicMock()
        mock_client.querymany.return_value = [
            {"query": "ENSG00000999999", "symbol": "LNC1"},
        ]

        with mock.patch.dict("sys.modules", {
            "mygene": mock.MagicMock(MyGeneInfo=mock.MagicMock(return_value=mock_client)),
        }):
            adata = _make_adata(
                ["ENSG00000000001", "ENSG00000999999", "ENSG00000999999"],
                var_columns={
                    "feature_name": ["TP53", "ENSG00000999999", "ENSG00000999999"],
                },
            )
            # 因为两个 ENSG00000999999 的 feature_name 相同，去重后会变成
            # ENSG00000999999 和 ENSG00000999999-1
            _sync_ensembl_to_symbol(adata)

        # 查询时传入的 ID 应已去掉 -1 后缀
        args, _kwargs = mock_client.querymany.call_args
        # 第一个 batch 的第一个参数是基因列表，应只含 ENSG00000999999（去后缀）
        query_batch = args[0]
        assert "ENSG00000999999-1" not in query_batch
        assert "ENSG00000999999" in query_batch


# ---------------------------------------------------------------------------
# _sync_ensembl_to_symbol — mygene 全量兜底（mock）
# ---------------------------------------------------------------------------


class TestSyncEnsemblToSymbolMygeneFallback:
    """无 feature_name 列时 mygene 全量兜底。"""

    def test_mygene_fallback_full(self):
        """无 feature_name → mygene 全量查询。"""
        mock_client = mock.MagicMock()
        mock_client.querymany.return_value = [
            {"query": "ENSG00000000001", "symbol": "TP53"},
            {"query": "ENSG00000000002", "symbol": "BRCA1"},
            {"query": "ENSG00000999999", "notfound": True},
        ]

        with mock.patch.dict("sys.modules", {
            "mygene": mock.MagicMock(MyGeneInfo=mock.MagicMock(return_value=mock_client)),
        }):
            adata = _make_adata(
                ["ENSG00000000001", "ENSG00000000002", "ENSG00000999999"],
            )
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                _sync_ensembl_to_symbol(adata)

        assert adata.var.index[0] == "TP53"
        assert adata.var.index[1] == "BRCA1"
        # notfound 的保留原 ENSG 作为索引
        assert adata.var.index[2] == "ENSG00000999999"
        assert adata.var["ensembl_id"].iloc[0] == "ENSG00000000001"
        assert adata.var["ensembl_id"].iloc[1] == "ENSG00000000002"
        # notfound 的基因：ensembl_id 清空（保留为索引但不重复存 ID）
        assert adata.var["ensembl_id"].iloc[2] == ""

    def test_mygene_fallback_dedup_after_query(self):
        """mygene 返回相同 symbol → 去重。"""
        mock_client = mock.MagicMock()
        mock_client.querymany.return_value = [
            {"query": "ENSG00000000001", "symbol": "SAME"},
            {"query": "ENSG00000000002", "symbol": "SAME"},
        ]

        with mock.patch.dict("sys.modules", {
            "mygene": mock.MagicMock(MyGeneInfo=mock.MagicMock(return_value=mock_client)),
        }):
            adata = _make_adata(
                ["ENSG00000000001", "ENSG00000000002"],
            )
            _sync_ensembl_to_symbol(adata)

        assert adata.var.index.is_unique
        assert adata.var.index[0] == "SAME"
        assert adata.var.index[1].startswith("SAME-")

    def test_mygene_fallback_network_error(self):
        """mygene 查询抛异常 → 保留 Ensembl ID，ensembl_id 留空。"""
        mock_client = mock.MagicMock()
        mock_client.querymany.side_effect = ConnectionError("offline")

        with mock.patch.dict("sys.modules", {
            "mygene": mock.MagicMock(MyGeneInfo=mock.MagicMock(return_value=mock_client)),
        }):
            adata = _make_adata(
                ["ENSG00000000001", "ENSG00000000002"],
            )
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                _sync_ensembl_to_symbol(adata)

        # 保留原 ENSG ID 作为索引
        assert adata.var.index[0] == "ENSG00000000001"
        assert adata.var.index[1] == "ENSG00000000002"
        # warning 已发出
        assert any("mygene 查询失败" in str(x.message) for x in w)


# ---------------------------------------------------------------------------
# _sync_ensembl_to_symbol — 幂等检查
# ---------------------------------------------------------------------------


class TestSyncEnsemblToSymbolIdempotent:
    def test_skips_when_ensembl_id_column_already_present(self):
        """ensembl_id 列已存在且非全空 → 幂等跳过。"""
        mock_mg_module = mock.MagicMock()
        mock_mg_module.MyGeneInfo = mock.MagicMock()

        adata = _make_adata(
            ["ENSG00000000001", "ENSG00000000002"],
            var_columns={
                "feature_name": ["TP53", "BRCA1"],
            },
        )
        # 预先填入 ensembl_id，模拟已处理过
        adata.var["ensembl_id"] = ["ENSG00000000001", "ENSG00000000002"]

        with mock.patch.dict("sys.modules", {"mygene": mock_mg_module}):
            _sync_ensembl_to_symbol(adata)
            # mygene 不应被调用（幂等检查直接返回）
            mock_mg_module.MyGeneInfo.assert_not_called()

        # 索引不变
        assert adata.var.index[0] == "ENSG00000000001"
        assert adata.var.index[1] == "ENSG00000000002"


# ---------------------------------------------------------------------------
# sync_gene_ids 集成：ensembl 输入 + feature_name
# ---------------------------------------------------------------------------


class TestSyncGeneIdsIntegration:
    def test_ensembl_index_with_feature_name_calls_sync_ensembl(self):
        """Ensembl 索引 + feature_name 列 → sync_gene_ids 走 ensembl 分支。"""
        mock_client = mock.MagicMock()
        # 二次补全查询全部 notfound，确保仍是 ENSG 的基因保留不变
        mock_client.querymany.return_value = [
            {"query": "ENSG00000000003", "notfound": True},
        ]

        with mock.patch.dict("sys.modules", {
            "mygene": mock.MagicMock(MyGeneInfo=mock.MagicMock(return_value=mock_client)),
        }):
            adata = _make_adata(
                ["ENSG00000000001", "ENSG00000000002", "ENSG00000000003"],
                var_columns={
                    "feature_name": ["TP53", "BRCA1", "ENSG00000000003"],
                },
            )
            sync_gene_ids(adata, "auto")

        assert adata.var.index[0] == "TP53"
        assert adata.var.index[1] == "BRCA1"
        assert adata.var.index[2] == "ENSG00000000003"
        assert adata.var["ensembl_id"].iloc[0] == "ENSG00000000001"

    def test_non_ensembl_index_does_not_trigger_sync_ensembl(self):
        """Symbol 索引 + gene_ids 列 → 不触发 _sync_ensembl_to_symbol，
        且 _sync_symbol_to_ensembl 直接使用 gene_ids，不调 mygene。"""
        mock_mg_module = mock.MagicMock()
        mock_mg_module.MyGeneInfo = mock.MagicMock()

        adata = _make_adata(
            ["TP53", "BRCA1", "MYC"],
            var_columns={"gene_ids": ["ENSG00000141510", "ENSG00000012048", "ENSG00000136997"]},
        )
        with mock.patch.dict("sys.modules", {"mygene": mock_mg_module}):
            sync_gene_ids(adata, "auto")
            # gene_ids 列有 Ensembl ID → _sync_symbol_to_ensembl 直接使用，不调 mygene
            mock_mg_module.MyGeneInfo.assert_not_called()

        assert "ensembl_id" in adata.var.columns
        assert adata.var["ensembl_id"].iloc[0] == "ENSG00000141510"
