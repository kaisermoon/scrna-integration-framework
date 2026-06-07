"""Tests for scrna_integration.io — read_with_manifest and helpers."""

from __future__ import annotations

import os
import tempfile
import warnings

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import yaml

from scrna_integration.io import (
    _compute_baseline_qc,
    _enforce_species,
    _generate_cell_id,
    _rename_original_annotations,
    _sync_gene_ids,
    _validate_manifest,
    _warn_layer2,
    read_with_manifest,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


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
        # Realistic ENSG IDs
        var_names = [f"ENSG{10000000000 + i}" for i in range(n_genes)]
    else:
        var_names = [f"GENE_{i}" for i in range(n_genes)]

    X = rng.poisson(5, size=(n_cells, n_genes)).astype(np.float32)
    if sparse:
        X = sp.csr_matrix(X)

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


class TestValidateManifest:
    def test_missing_species_raises(self):
        m = _make_manifest_yaml()
        del m["species"]
        with pytest.raises(ValueError, match="species"):
            _validate_manifest(m)

    def test_non_human_species_raises(self):
        m = _make_manifest_yaml({"species": "mouse"})
        with pytest.raises(ValueError, match="species"):
            _validate_manifest(m)

    def test_missing_input_raises(self):
        m = _make_manifest_yaml()
        del m["input"]
        with pytest.raises(ValueError, match="input"):
            _validate_manifest(m)

    def test_missing_input_format_raises(self):
        m = _make_manifest_yaml({"input": {"path": "/tmp/test.h5ad"}})
        with pytest.raises(ValueError, match="format"):
            _validate_manifest(m)

    def test_invalid_input_format_raises(self):
        m = _make_manifest_yaml({"input": {"format": "excel", "path": "/tmp/x.xlsx"}})
        with pytest.raises(ValueError, match="Unsupported"):
            _validate_manifest(m)

    def test_missing_source_dataset_raises(self):
        m = _make_manifest_yaml()
        del m["source_dataset"]
        with pytest.raises(ValueError, match="source_dataset"):
            _validate_manifest(m)

    def test_missing_project_id_raises(self):
        m = _make_manifest_yaml()
        del m["project_id"]
        with pytest.raises(ValueError, match="project_id"):
            _validate_manifest(m)

    def test_missing_disease_system_raises(self):
        m = _make_manifest_yaml()
        del m["disease_system"]
        with pytest.raises(ValueError, match="disease_system"):
            _validate_manifest(m)

    def test_missing_original_annotations_raises(self):
        m = _make_manifest_yaml()
        del m["original_annotations"]
        with pytest.raises(ValueError, match="original_annotations"):
            _validate_manifest(m)

    def test_qc_overrides_skip_without_reason_raises(self):
        m = _make_manifest_yaml({"qc_overrides": {"doublet_removal": {"skip": True}}})
        with pytest.raises(ValueError, match="reason"):
            _validate_manifest(m)

    def test_qc_overrides_skip_with_reason_passes(self):
        m = _make_manifest_yaml(
            {"qc_overrides": {"doublet_removal": {"skip": True, "reason": "too few cells"}}}
        )
        _validate_manifest(m)  # should not raise

    def test_valid_minimal_manifest_passes(self):
        m = _make_manifest_yaml()
        _validate_manifest(m)


# ---------------------------------------------------------------------------
# Species enforcement
# ---------------------------------------------------------------------------


class TestEnforceSpecies:
    def test_human_accepted(self):
        adata = _make_synthetic_adata()
        _enforce_species(adata, "human")
        assert adata.uns["species"] == "human"

    def test_non_human_raises(self):
        adata = _make_synthetic_adata()
        with pytest.raises(ValueError, match="species.*not supported"):
            _enforce_species(adata, "mouse")


# ---------------------------------------------------------------------------
# Cell ID generation
# ---------------------------------------------------------------------------


class TestGenerateCellId:
    def test_format(self):
        adata = _make_synthetic_adata(n_cells=10, n_genes=5)
        _generate_cell_id(adata, "TestDS")
        assert all(adata.obs["cell_id"].str.startswith("TestDS_"))
        assert adata.obs["cell_id"].iloc[0] == "TestDS_sample_0_cell_0"

    def test_unique_cell_ids(self):
        adata = _make_synthetic_adata(n_cells=100)
        _generate_cell_id(adata, "DS")
        assert adata.obs["cell_id"].nunique() == 100


# ---------------------------------------------------------------------------
# Original annotations rename
# ---------------------------------------------------------------------------


class TestRenameOriginalAnnotations:
    def test_rename_with_role(self):
        adata = _make_synthetic_adata()
        adata.obs["author_label"] = ["T_cell", "B_cell"] * 50
        manifest = {
            "source_dataset": "Test_2024",
            "original_annotations": [
                {"column": "author_label", "role": "primary", "granularity": "broad"}
            ],
        }
        _rename_original_annotations(adata, manifest)
        assert "cell_type_original_Test_2024_v1_primary" in adata.obs.columns
        assert adata.obs["cell_type_original_Test_2024_v1_primary"].iloc[0] == "T_cell"

    def test_rename_without_role(self):
        adata = _make_synthetic_adata()
        adata.obs["label"] = ["type_A"] * 100
        manifest = {
            "source_dataset": "DS",
            "original_annotations": [{"column": "label"}],
        }
        _rename_original_annotations(adata, manifest)
        assert "cell_type_original_DS_v1" in adata.obs.columns

    def test_missing_column_warns(self):
        adata = _make_synthetic_adata()
        manifest = {
            "source_dataset": "DS",
            "original_annotations": [{"column": "nonexistent"}],
        }
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _rename_original_annotations(adata, manifest)
            assert len(w) == 1
            assert "nonexistent" in str(w[0].message)


# ---------------------------------------------------------------------------
# Gene ID sync
# ---------------------------------------------------------------------------


class TestSyncGeneIds:
    def test_symbol_index_with_gene_ids_column(self):
        """When var has gene_ids column with Ensembl IDs, use it directly."""
        adata = _make_synthetic_adata(n_genes=10, gene_index="symbol")
        # Add 10x-style gene_ids column
        adata.var["gene_ids"] = [f"ENSG{20000000000 + i}" for i in range(10)]
        _sync_gene_ids(adata, "auto")
        assert "ensembl_id" in adata.var.columns
        assert adata.var["ensembl_id"].iloc[0] == "ENSG20000000000"
        # var.index should remain symbols
        assert adata.var.index[0] == "GENE_0"

    def test_symbol_index_no_gene_ids_falls_back_to_mygene(self):
        """When no gene_ids column, try mygene.  At minimum, ensembl_id column is created."""
        adata = _make_synthetic_adata(n_genes=5, gene_index="symbol")
        # Use synthetic gene names that won't be found by mygene
        _sync_gene_ids(adata, "auto")
        assert "ensembl_id" in adata.var.columns
        # Column should exist (may be empty for synthetic genes)
        assert adata.var["ensembl_id"].isna().sum() >= 0

    def test_ensembl_index_with_feature_name(self):
        """Ensembl var.index with feature_name column -> use feature_name as index."""
        # Start with symbols as index then replace with ensembl
        adata = _make_synthetic_adata(n_genes=5, gene_index="symbol")
        adata.var["feature_name"] = adata.var.index.values  # Original symbols
        adata.var.index = [f"ENSG{30000000000 + i}" for i in range(5)]
        _sync_gene_ids(adata, "auto")
        assert not adata.var.index[0].startswith("ENSG")
        assert "ensembl_id" in adata.var.columns

    def test_ensembl_index_via_mygene(self):
        """Ensembl index without feature_name -> try mygene."""
        adata = _make_synthetic_adata(n_genes=3, gene_index="ensembl")
        # Use like-looking but synthetic ensembl IDs that mygene won't find
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _sync_gene_ids(adata, "auto")
        assert "ensembl_id" in adata.var.columns


# ---------------------------------------------------------------------------
# Baseline QC
# ---------------------------------------------------------------------------


class TestBaselineQC:
    def test_adds_qc_columns(self):
        adata = _make_synthetic_adata(n_cells=20, n_genes=50)
        # Add a couple MT genes for testing
        new_var = list(adata.var.index)
        new_var[0] = "MT-CO1"
        new_var[5] = "MT-ND1"
        new_var[10] = "RPS3"
        adata.var.index = new_var
        _compute_baseline_qc(adata)
        assert "n_genes" in adata.obs.columns
        assert "total_counts" in adata.obs.columns
        assert "pct_counts_mt" in adata.obs.columns
        assert "pct_counts_ribo" in adata.obs.columns
        assert adata.obs["n_genes"].min() >= 0
        assert adata.obs["total_counts"].min() >= 0

    def test_no_mt_genes_no_pct_mt(self):
        adata = _make_synthetic_adata(n_cells=10, n_genes=10)
        _compute_baseline_qc(adata)
        assert "pct_counts_mt" not in adata.obs.columns


# ---------------------------------------------------------------------------
# Layer 2 warn
# ---------------------------------------------------------------------------


class TestWarnLayer2:
    def test_missing_fields_emit_warnings(self):
        adata = _make_synthetic_adata()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _warn_layer2(adata)
            # Should warn for all 7 missing fields
            assert len(w) == 7
            for warning in w:
                assert "Layer2" in str(warning.message)

    def test_does_not_raise(self):
        """Layer 2 warnings should never raise or block."""
        adata = _make_synthetic_adata()
        _warn_layer2(adata)  # should not raise


# ---------------------------------------------------------------------------
# RDS format
# ---------------------------------------------------------------------------


class TestRDS:
    def test_not_implemented(self):
        manifest = _make_manifest_yaml(
            {"input": {"format": "rds", "path": "/tmp/test.rds"}}
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = _write_manifest(manifest, tmpdir)
            with pytest.raises(NotImplementedError, match="RDS"):
                read_with_manifest(path)


# ---------------------------------------------------------------------------
# End-to-end: synthetic h5ad
# ---------------------------------------------------------------------------


class TestReadWithManifestSynthetic:
    def test_minimal_manifest_reads_h5ad(self):
        """Full end-to-end with a synthetic h5ad."""
        adata = _make_synthetic_adata()
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)

            manifest = _make_manifest_yaml(
                {"input": {"format": "h5ad", "path": h5ad_path}}
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)

        assert isinstance(result, anndata.AnnData)
        assert result.uns["species"] == "human"
        assert "disease_system" in result.obs.columns
        assert "cell_id" in result.obs.columns
        assert "ensembl_id" in result.var.columns
        assert "n_genes" in result.obs.columns
        assert result.obs["disease_system"].iloc[0] == "gastric"

    def test_obs_mapping_applied(self):
        """obs_mapping renames columns and value_mapping transforms values."""
        adata = _make_synthetic_adata(n_cells=10)
        adata.obs["orig_status"] = ["CAG", "IM", "normal"] * 3 + ["CAG"]
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {
                    "input": {"format": "h5ad", "path": h5ad_path},
                    "obs_mapping": {"disease": "orig_status"},
                    "value_mapping": {
                        "disease": {"CAG": "atrophic_gastritis", "IM": "metaplasia"}
                    },
                }
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)

        assert "disease" in result.obs.columns
        # CAG -> atrophic_gastritis, IM -> metaplasia, normal unchanged
        values = set(result.obs["disease"])
        assert "atrophic_gastritis" in values
        assert "metaplasia" in values

    def test_original_annotations_renamed(self):
        """Author annotation columns get renamed with source_dataset prefix."""
        adata = _make_synthetic_adata()
        adata.obs["CellType"] = ["T", "B"] * 50
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {
                    "input": {"format": "h5ad", "path": h5ad_path},
                    "original_annotations": [
                        {"column": "CellType", "role": "primary"}
                    ],
                }
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)

        assert "cell_type_original_Test_2024_v1_primary" in result.obs.columns

    def test_species_non_human_fails(self):
        """Non-human species should raise ValueError."""
        adata = _make_synthetic_adata()
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {"input": {"format": "h5ad", "path": h5ad_path}, "species": "mouse"}
            )
            mpath = _write_manifest(manifest, tmpdir)
            with pytest.raises(ValueError, match="species"):
                read_with_manifest(mpath)

    def test_raw_path_recorded(self):
        """raw_path in input block -> adata.uns['raw_matrix_path']."""
        adata = _make_synthetic_adata()
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {
                    "input": {
                        "format": "h5ad",
                        "path": h5ad_path,
                        "raw_path": "/some/raw/path",
                    }
                }
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)
        assert result.uns["raw_matrix_path"] == "/some/raw/path"

    def test_no_raw_path_gives_none(self):
        adata = _make_synthetic_adata()
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {"input": {"format": "h5ad", "path": h5ad_path}}
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)
        assert result.uns["raw_matrix_path"] is None


# ---------------------------------------------------------------------------
# 10x_mtx synthetic test
# ---------------------------------------------------------------------------


class TestRead10xMtx:
    def test_reads_single_10x_dir(self):
        """Read a 10x mtx directory written by scanpy."""
        adata = _make_synthetic_adata(n_cells=30, n_genes=10)
        with tempfile.TemporaryDirectory() as tmpdir:
            mtx_dir = os.path.join(tmpdir, "filtered_feature_bc_matrix")
            os.makedirs(mtx_dir)
            # Write in 10x format (scanpy expects .mtx.gz)
            import gzip
            import io as sys_io

            scipy_io = pytest.importorskip("scipy.io")
            buf = sys_io.BytesIO()
            scipy_io.mmwrite(buf, adata.X.T)  # 10x convention: genes x cells
            with gzip.open(os.path.join(mtx_dir, "matrix.mtx.gz"), "wb") as f:
                f.write(buf.getvalue())
            # features.tsv.gz
            with gzip.open(os.path.join(mtx_dir, "features.tsv.gz"), "wt") as f:
                for g in adata.var.index:
                    f.write(f"{g}\t{g}\tGene Expression\n")
            # barcodes.tsv.gz
            with gzip.open(os.path.join(mtx_dir, "barcodes.tsv.gz"), "wt") as f:
                for bc in adata.obs_names:
                    f.write(f"{bc}\n")

            manifest = _make_manifest_yaml(
                {"input": {"format": "10x_mtx", "path": mtx_dir}}
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)

        assert isinstance(result, anndata.AnnData)
        assert result.shape[0] == 30
        assert "sample_id" in result.obs.columns


# ---------------------------------------------------------------------------
# txt.gz synthetic test
# ---------------------------------------------------------------------------


class TestReadTxtGz:
    def test_reads_txt_gz_directory(self):
        """Read tab-separated gzipped count matrix."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a txt.gz file in gene x cells layout
            import gzip

            cells = ["cell_1", "cell_2", "cell_3"]
            fp = os.path.join(tmpdir, "sample1_count.txt.gz")
            with gzip.open(fp, "wt") as f:
                f.write("gene\t" + "\t".join(cells) + "\n")
                f.write("GENE_A\t1\t2\t0\n")
                f.write("GENE_B\t0\t3\t1\n")
                f.write("GENE_C\t4\t0\t2\n")

            manifest = _make_manifest_yaml(
                {"input": {"format": "txt.gz", "path": tmpdir}}
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)

        assert result.shape == (3, 3)
        assert result.obs["sample_id"].iloc[0] == "sample1"


# ---------------------------------------------------------------------------
# H5 format (10x .h5)
# ---------------------------------------------------------------------------


class TestReadH5:
    def test_reads_10x_h5(self):
        """Read a 10x .h5 file — requires real 10x HDF5 format fixture.

        Synthetic h5ad files cannot substitute for the 10x HDF5 format
        (different HDF5 structure).  The end-to-end tests on real datasets
        cover the h5 path when a 10x h5 fixture is available.
        """
        pytest.skip("10x H5 format requires a real 10x h5 fixture (not h5ad)")
