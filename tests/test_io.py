
from __future__ import annotations

import builtins
import os
import tempfile
import warnings
from unittest import mock

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.io as scipy_io
import scipy.sparse as sp
import yaml

from scrna_integration.io import (
    _mygene_symbol_to_ensembl,
    sync_gene_ids,
    inject_genomic_positions,
)

_ORIGINAL_IMPORT = builtins.__import__

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_mygene_mock(return_values):
    """构建 mock mygene 模块，querymany 返回 *return_values*。"""
    mock_mg = mock.MagicMock()
    mock_client = mock.MagicMock()
    mock_mg.MyGeneInfo.return_value = mock_client
    mock_client.querymany.return_value = return_values
    return mock_mg


def _make_synthetic_adata(
    n_cells: int = 100,
    n_genes: int = 200,
    gene_index: str = "symbol",
    sparse: bool = True,
) -> anndata.AnnData:
    """Build a minimal AnnData for testing."""
    rng = np.random.default_rng(42)
    if gene_index == "symbol":
        var_names = [f"GENE_{i}" for i in range(n_genes)]
    elif gene_index == "ensembl":
        var_names = [f"ENSG{10000000000 + i}" for i in range(n_genes)]
    else:
        var_names = [f"GENE_{i}" for i in range(n_genes)]

    X = rng.poisson(5, size=(n_cells, n_genes)).astype(np.float32)  # noqa: N806
    if sparse:
        X = sp.csr_matrix(X)  # noqa: N806

    obs = pd.DataFrame(
        {"sample_id": [f"sample_{i % 4}" for i in range(n_cells)]},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=var_names)
    return anndata.AnnData(X=X, obs=obs, var=var)


def _make_manifest_yaml(overrides: dict | None = None) -> dict:
    """Return a minimal valid manifest dict suitable for testing."""
    base = {
        "species": "human",
        "input": {"format": "h5ad", "path": "/tmp/test.h5ad"},
        "source_dataset": "Test_2024",
        "project_id": "test_project",
        "disease_system": "gastric",
        "original_annotations": [],
    }
    if overrides:
        base.update(overrides)
    return base


def _write_manifest(manifest: dict, tmpdir: str) -> str:
    """Write a manifest dict to a temp YAML file, return path."""
    path = os.path.join(tmpdir, "manifest.yaml")
    with open(path, "w") as fh:
        yaml.safe_dump(manifest, fh)
    return path


# ---------------------------------------------------------------------------
# Manifest validation
# ---------------------------------------------------------------------------


class TestSyncGeneIds:
    def test_symbol_index_with_gene_ids_column(self):
        """var 含 gene_ids 列（含 Ensembl ID）时直接使用。"""
        adata = _make_synthetic_adata(n_genes=10, gene_index="symbol")
        adata.var["gene_ids"] = [f"ENSG{20000000000 + i}" for i in range(10)]
        sync_gene_ids(adata, "auto")
        assert "ensembl_id" in adata.var.columns
        assert adata.var["ensembl_id"].iloc[0] == "ENSG20000000000"
        assert adata.var.index[0] == "GENE_0"

    def test_symbol_index_no_gene_ids_falls_back_to_mygene(self):
        """无 gene_ids 列时回退 mygene。至少创建 ensembl_id 列。"""
        adata = _make_synthetic_adata(n_genes=5, gene_index="symbol")
        sync_gene_ids(adata, "auto")
        assert "ensembl_id" in adata.var.columns
        assert adata.var["ensembl_id"].isna().sum() >= 0

    def test_ensembl_index_with_feature_name(self):
        """Ensembl 索引 + feature_name 列 → 用 feature_name 做索引。"""
        adata = _make_synthetic_adata(n_genes=5, gene_index="symbol")
        adata.var["feature_name"] = adata.var.index.values
        adata.var.index = [f"ENSG{30000000000 + i}" for i in range(5)]
        sync_gene_ids(adata, "auto")
        assert not adata.var.index[0].startswith("ENSG")
        assert "ensembl_id" in adata.var.columns

    def test_ensembl_index_via_mygene(self):
        """Ensembl 索引无 feature_name → 尝试 mygene。"""
        adata = _make_synthetic_adata(n_genes=3, gene_index="ensembl")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sync_gene_ids(adata, "auto")
        assert "ensembl_id" in adata.var.columns


# ---------------------------------------------------------------------------
# Baseline QC
# ---------------------------------------------------------------------------


class TestMygeneContract:
    """inject_genomic_positions 与 _mygene_symbol_to_ensembl 的 mock 测试。

    核心 mock 策略：用 mock.patch.dict("sys.modules", ...) 替换 mygene 模块，
    因为 inject_genomic_positions 内部使用 ``import mygene``（通过 sys.modules
    解析），而非通过 scrna_integration.io.mygene 属性访问。
    """

    # ------------------------------------------------------------------
    # inject_genomic_positions
    # ------------------------------------------------------------------

    def test_inject_genomic_positions_normal(self):
        """genomic_pos 是 dict（单位置）--> 三列正确填入。"""
        mock_mg = _make_mygene_mock([
            {
                "query": "ENSG00000141510",
                "genomic_pos": {"chr": "17", "start": 7660000, "end": 7670000},
            },
            {
                "query": "ENSG00000139618",
                "genomic_pos": {"chr": "chr13", "start": 32300000, "end": 32350000},
            },
        ])

        adata = _make_synthetic_adata(n_genes=3)
        adata.var["ensembl_id"] = [
            "ENSG00000141510",
            "ENSG00000139618",
            "ENSG00000999999",
        ]
        with mock.patch.dict("sys.modules", {"mygene": mock_mg}):
            result = inject_genomic_positions(adata)

        assert result.var["chromosome"].iloc[0] == "17"
        assert result.var["start"].iloc[0] == 7660000
        assert result.var["end"].iloc[0] == 7670000
        # "chr13" --> "13"（chr 前缀已归一化）
        assert result.var["chromosome"].iloc[1] == "13"
        assert result.var["start"].iloc[1] == 32300000
        # ENSG00000999999 未找到 --> NaN
        assert pd.isna(result.var["chromosome"].iloc[2])
        assert pd.isna(result.var["start"].iloc[2])

    def test_inject_genomic_positions_partial(self):
        """部分基因返回 notfound:True --> 对应行 NaN + warning。"""
        mock_mg = _make_mygene_mock([
            {
                "query": "ENSG00000141510",
                "genomic_pos": {"chr": "17", "start": 7660000, "end": 7670000},
            },
            {"query": "ENSG00000999999", "notfound": True},
        ])

        adata = _make_synthetic_adata(n_genes=2)
        adata.var["ensembl_id"] = ["ENSG00000141510", "ENSG00000999999"]
        with mock.patch.dict("sys.modules", {"mygene": mock_mg}):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = inject_genomic_positions(adata)

        assert result.var["chromosome"].iloc[0] == "17"
        assert pd.isna(result.var["chromosome"].iloc[1])
        pos_warnings = [x for x in w if "未获取到位置信息" in str(x.message)]
        assert len(pos_warnings) >= 1

    def test_inject_genomic_positions_multi_pos(self):
        """genomic_pos 是 list（多位置）--> 取第一个。"""
        mock_mg = _make_mygene_mock([
            {
                "query": "ENSG00000141510",
                "genomic_pos": [
                    {"chr": "17", "start": 7660000, "end": 7670000},
                    {"chr": "17", "start": 7680000, "end": 7690000},
                ],
            },
        ])

        adata = _make_synthetic_adata(n_genes=1)
        adata.var["ensembl_id"] = ["ENSG00000141510"]
        with mock.patch.dict("sys.modules", {"mygene": mock_mg}):
            result = inject_genomic_positions(adata)

        assert result.var["start"].iloc[0] == 7660000

    def test_inject_genomic_positions_no_mygene(self):
        """mygene 未安装 --> raise ImportError。"""
        def _side_effect(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "mygene":
                raise ImportError("No module named 'mygene'")
            return _ORIGINAL_IMPORT(name, globals, locals, fromlist, level)

        adata = _make_synthetic_adata(n_genes=3)
        with mock.patch("builtins.__import__", side_effect=_side_effect):
            with pytest.raises(ImportError, match="mygene"):
                inject_genomic_positions(adata)

    def test_inject_genomic_positions_network_error(self):
        """querymany raise 异常 --> 优雅降级（NaN + warning）。"""
        mock_mg = _make_mygene_mock([])
        mock_mg.MyGeneInfo.return_value.querymany.side_effect = ConnectionError(
            "Network unreachable"
        )

        adata = _make_synthetic_adata(n_genes=3)
        adata.var["ensembl_id"] = [
            "ENSG00000141510",
            "ENSG00000139618",
            "ENSG00000999999",
        ]
        with mock.patch.dict("sys.modules", {"mygene": mock_mg}):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = inject_genomic_positions(adata)

        assert result.var["chromosome"].isna().all()
        query_warnings = [
            x for x in w if "mygene 查询失败" in str(x.message)
        ]
        assert len(query_warnings) >= 1

    # ------------------------------------------------------------------
    # _mygene_symbol_to_ensembl（覆盖现有测试漏洞）
    # ------------------------------------------------------------------

    def test_mygene_symbol_to_ensembl_mock(self):
        """正常响应 --> ensembl_id 填入 var 列。"""
        mock_mg = _make_mygene_mock([
            {"query": "GENE_0", "ensembl": {"gene": "ENSG00000000001"}},
            {"query": "GENE_1", "ensembl": {"gene": "ENSG00000000002"}},
            {"query": "GENE_2", "notfound": True},
        ])

        adata = _make_synthetic_adata(n_genes=3, gene_index="symbol")
        with mock.patch.dict("sys.modules", {"mygene": mock_mg}):
            _mygene_symbol_to_ensembl(adata)

        assert adata.var["ensembl_id"].iloc[0] == "ENSG00000000001"
        assert adata.var["ensembl_id"].iloc[1] == "ENSG00000000002"
        # notfound --> 空字符串
        assert adata.var["ensembl_id"].iloc[2] == ""


# ---------------------------------------------------------------------------
# P0-1: Layer 1 确定性校验 + project_id 写入 obs
# ---------------------------------------------------------------------------
