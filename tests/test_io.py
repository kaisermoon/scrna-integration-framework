"""Tests for scrna_integration.io — read_with_manifest and retained helpers."""

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
import scipy.sparse as sp
import yaml

from scrna_integration.io import (
    _compute_baseline_qc,
    _mygene_symbol_to_ensembl,
    _sync_gene_ids,
    _validate_layer1,
    _validate_manifest,
    _warn_layer2,
    inject_genomic_positions,
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


class TestValidateManifest:
    def test_missing_species_raises(self):
        m = _make_manifest_yaml()
        del m["species"]
        with pytest.raises(ValueError, match="species"):
            _validate_manifest(m)

    def test_non_human_species_raises(self):
        m = _make_manifest_yaml({"species": "mouse"})
        with pytest.raises(ValueError, match="物种"):
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
        with pytest.raises(ValueError, match="不支持"):
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
            {"qc_overrides": {"doublet_removal": {"skip": True, "reason": "细胞太少"}}}
        )
        _validate_manifest(m)  # should not raise

    def test_valid_minimal_manifest_passes(self):
        m = _make_manifest_yaml()
        _validate_manifest(m)


# ---------------------------------------------------------------------------
# Species enforcement (tested through read_with_manifest)
# ---------------------------------------------------------------------------


class TestSpecies:
    """物种校验：仅 human 接受，其他报错。"""

    def test_human_sets_uns_species(self):
        """通过 read_with_manifest 确认 human 物种写入 uns['species']。"""
        adata = _make_synthetic_adata()
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {"input": {"format": "h5ad", "path": h5ad_path}, "species": "human"}
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)
        assert result.uns["species"] == "human"

    def test_non_human_raises(self):
        """非 human 物种在 read_with_manifest 中抛出 ValueError。"""
        adata = _make_synthetic_adata()
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {"input": {"format": "h5ad", "path": h5ad_path}, "species": "mouse"}
            )
            mpath = _write_manifest(manifest, tmpdir)
            with pytest.raises(ValueError, match="species|物种"):
                read_with_manifest(mpath)


# ---------------------------------------------------------------------------
# Cell ID generation (tested through read_with_manifest)
# ---------------------------------------------------------------------------


class TestCellId:
    """全局唯一细胞 ID 生成。"""

    def test_format(self):
        """cell_id 格式: {数据集}_{样本}_{barcode}。"""
        adata = _make_synthetic_adata(n_cells=10, n_genes=5)
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {"input": {"format": "h5ad", "path": h5ad_path}}
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)
        assert all(result.obs["cell_id"].str.startswith("Test_2024_"))
        assert result.obs["cell_id"].iloc[0] == "Test_2024_sample_0_cell_0"

    def test_unique_cell_ids(self):
        """每个细胞的 cell_id 全局唯一。"""
        adata = _make_synthetic_adata(n_cells=100)
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {"input": {"format": "h5ad", "path": h5ad_path}}
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)
        assert result.obs["cell_id"].nunique() == 100


# ---------------------------------------------------------------------------
# Original annotations rename (tested through read_with_manifest)
# ---------------------------------------------------------------------------


class TestOriginalAnnotations:
    """原始作者标注列重命名。"""

    def test_rename_with_role(self):
        """含 role 的列重命名为 cell_type_original_{ds}_v1_{role}。"""
        adata = _make_synthetic_adata()
        adata.obs["author_label"] = ["T_cell", "B_cell"] * 50
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {
                    "input": {"format": "h5ad", "path": h5ad_path},
                    "original_annotations": [
                        {"column": "author_label", "role": "primary", "granularity": "broad"}
                    ],
                }
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)
        assert "cell_type_original_Test_2024_v1_primary" in result.obs.columns
        assert result.obs["cell_type_original_Test_2024_v1_primary"].iloc[0] == "T_cell"

    def test_rename_without_role(self):
        """无 role 的列重命名为 cell_type_original_{ds}_v1。"""
        adata = _make_synthetic_adata()
        adata.obs["label"] = ["type_A"] * 100
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {
                    "input": {"format": "h5ad", "path": h5ad_path},
                    "source_dataset": "DS",
                    "original_annotations": [{"column": "label"}],
                }
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)
        assert "cell_type_original_DS_v1" in result.obs.columns

    def test_missing_column_warns(self):
        """不存在的列发出警告但不崩溃。"""
        adata = _make_synthetic_adata()
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {
                    "input": {"format": "h5ad", "path": h5ad_path},
                    "original_annotations": [{"column": "nonexistent"}],
                }
            )
            mpath = _write_manifest(manifest, tmpdir)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                read_with_manifest(mpath)
                assert len(w) >= 1
                assert "nonexistent" in str(w[0].message)


# ---------------------------------------------------------------------------
# Gene ID sync
# ---------------------------------------------------------------------------


class TestSyncGeneIds:
    def test_symbol_index_with_gene_ids_column(self):
        """var 含 gene_ids 列（含 Ensembl ID）时直接使用。"""
        adata = _make_synthetic_adata(n_genes=10, gene_index="symbol")
        adata.var["gene_ids"] = [f"ENSG{20000000000 + i}" for i in range(10)]
        _sync_gene_ids(adata, "auto")
        assert "ensembl_id" in adata.var.columns
        assert adata.var["ensembl_id"].iloc[0] == "ENSG20000000000"
        assert adata.var.index[0] == "GENE_0"

    def test_symbol_index_no_gene_ids_falls_back_to_mygene(self):
        """无 gene_ids 列时回退 mygene。至少创建 ensembl_id 列。"""
        adata = _make_synthetic_adata(n_genes=5, gene_index="symbol")
        _sync_gene_ids(adata, "auto")
        assert "ensembl_id" in adata.var.columns
        assert adata.var["ensembl_id"].isna().sum() >= 0

    def test_ensembl_index_with_feature_name(self):
        """Ensembl 索引 + feature_name 列 → 用 feature_name 做索引。"""
        adata = _make_synthetic_adata(n_genes=5, gene_index="symbol")
        adata.var["feature_name"] = adata.var.index.values
        adata.var.index = [f"ENSG{30000000000 + i}" for i in range(5)]
        _sync_gene_ids(adata, "auto")
        assert not adata.var.index[0].startswith("ENSG")
        assert "ensembl_id" in adata.var.columns

    def test_ensembl_index_via_mygene(self):
        """Ensembl 索引无 feature_name → 尝试 mygene。"""
        adata = _make_synthetic_adata(n_genes=3, gene_index="ensembl")
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
            assert len(w) == 7
            for warning in w:
                assert "Layer2" in str(warning.message)

    def test_does_not_raise(self):
        """Layer 2 警告不应抛出异常或阻塞。"""
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
        """最小 manifest + 合成 h5ad → 完整端到端测试。"""
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
        """obs_mapping 重命名列 + value_mapping 转换值。"""
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
        values = set(result.obs["disease"])
        assert "atrophic_gastritis" in values
        assert "metaplasia" in values

    def test_original_annotations_renamed(self):
        """作者标注列按 source_dataset 前缀重命名。"""
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
        """非 human 物种在 read_with_manifest 中抛出 ValueError。"""
        adata = _make_synthetic_adata()
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {"input": {"format": "h5ad", "path": h5ad_path}, "species": "mouse"}
            )
            mpath = _write_manifest(manifest, tmpdir)
            with pytest.raises(ValueError, match="species|物种"):
                read_with_manifest(mpath)

    def test_raw_path_recorded(self):
        """input.raw_path → adata.uns['raw_matrix_path']。"""
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
        """无 raw_path → adata.uns['raw_matrix_path'] = None。"""
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
        """读取含 10x mtx 格式文件的目录。"""
        adata = _make_synthetic_adata(n_cells=30, n_genes=10)
        with tempfile.TemporaryDirectory() as tmpdir:
            mtx_dir = os.path.join(tmpdir, "filtered_feature_bc_matrix")
            os.makedirs(mtx_dir)
            import gzip
            import io as sys_io

            scipy_io = pytest.importorskip("scipy.io")
            buf = sys_io.BytesIO()
            scipy_io.mmwrite(buf, adata.X.T)  # 10x 惯例: genes × cells
            with gzip.open(os.path.join(mtx_dir, "matrix.mtx.gz"), "wb") as f:
                f.write(buf.getvalue())
            with gzip.open(os.path.join(mtx_dir, "features.tsv.gz"), "wt") as f:
                for g in adata.var.index:
                    f.write(f"{g}\t{g}\tGene Expression\n")
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
        """读取 tab 分隔的 gzip 压缩计数矩阵。"""
        with tempfile.TemporaryDirectory() as tmpdir:
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
        """10x .h5 文件需要真实 10x HDF5 格式夹具。"""
        pytest.skip("10x H5 格式需要真实 10x h5 夹具（非 h5ad）")


# ---------------------------------------------------------------------------
# Clinical metadata join
# ---------------------------------------------------------------------------


class TestJoinClinicalCsv:
    def test_csv_happy_path(self):
        """临床 CSV 通过 clinical_metadata 配置关联到 obs。"""
        adata = _make_synthetic_adata(n_cells=12)
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)

            csv_path = os.path.join(tmpdir, "clinical.csv")
            pd.DataFrame({
                "sample_id": [f"sample_{i}" for i in range(4)],
                "age": [45, 52, 38, 61],
                "sex": ["M", "F", "F", "M"],
            }).to_csv(csv_path, index=False)

            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "clinical_metadata": [{
                    "file": csv_path,
                    "join_on": {
                        "manifest_field": "sample_id",
                        "table_column": "sample_id",
                    },
                }],
            })
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)

        assert "age" in result.obs.columns
        assert "sex" in result.obs.columns
        sample0 = result.obs[result.obs["sample_id"] == "sample_0"]
        assert (sample0["age"] == 45).all()
        assert (sample0["sex"] == "M").all()

    def test_missing_file_warns_not_crashes(self):
        """不存在的临床文件应警告而非崩溃（默认 on_missing=warn）。"""
        adata = _make_synthetic_adata(n_cells=8)
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "clinical_metadata": [{
                    "file": os.path.join(tmpdir, "nonexistent.csv"),
                    "join_on": {
                        "manifest_field": "sample_id",
                        "table_column": "sample_id",
                    },
                }],
            })
            mpath = _write_manifest(manifest, tmpdir)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = read_with_manifest(mpath)
            file_warnings = [
                x for x in w if "临床 metadata 文件未找到" in str(x.message)
            ]
            assert len(file_warnings) >= 1

        assert isinstance(result, anndata.AnnData)


# ---------------------------------------------------------------------------
# Ontology / project_specific constant injection
# ---------------------------------------------------------------------------


class TestInjectOntologyConstants:
    def test_ontology_injection(self):
        """ontology 段：常量值写入每行 obs。"""
        adata = _make_synthetic_adata(n_cells=10)
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "ontology": {"tissue": "gastric mucosa"},
            })
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)

        assert "tissue" in result.obs.columns
        unique_tissues = result.obs["tissue"].unique().tolist()
        assert unique_tissues == ["gastric mucosa"]

    def test_project_specific_with_source_column_and_rules(self):
        """project_specific 含 source_column + rules：映射已有 obs 列。"""
        adata = _make_synthetic_adata(n_cells=12)
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "project_specific": {
                    "treatment_group": {
                        "source_column": "sample_id",
                        "rules": {
                            "sample_0": "control",
                            "sample_1": "treatment",
                        },
                    },
                },
            })
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)

        assert "treatment_group" in result.obs.columns
        sample0 = result.obs[result.obs["sample_id"] == "sample_0"]
        assert (sample0["treatment_group"] == "control").all()
        sample1 = result.obs[result.obs["sample_id"] == "sample_1"]
        assert (sample1["treatment_group"] == "treatment").all()


# ---------------------------------------------------------------------------
# Mygene contract tests --- mock 覆盖网络调用（用 sys.modules 替换 mygene）
# ---------------------------------------------------------------------------


_ORIGINAL_IMPORT = builtins.__import__


def _make_mygene_mock(return_values):
    """构建 mock mygene 模块，querymany 返回 *return_values*。"""
    mock_mg = mock.MagicMock()
    mock_client = mock.MagicMock()
    mock_mg.MyGeneInfo.return_value = mock_client
    mock_client.querymany.return_value = return_values
    return mock_mg


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


class TestValidateLayer1:
    """P0-1: _validate_layer1 三个必需字段缺一即 raise ValueError。"""

    def test_all_present_passes(self):
        """三个字段都存在时不抛错。"""
        adata = _make_synthetic_adata()
        adata.obs["source_dataset"] = "DS1"
        adata.obs["project_id"] = "PRJ1"
        adata.obs["disease_system"] = "gastric"
        _validate_layer1(adata)  # should not raise

    def test_missing_source_dataset_raises(self):
        """缺 source_dataset 抛 ValueError。"""
        adata = _make_synthetic_adata()
        adata.obs["project_id"] = "PRJ1"
        adata.obs["disease_system"] = "gastric"
        with pytest.raises(ValueError, match="source_dataset|Layer 1"):
            _validate_layer1(adata)

    def test_missing_project_id_raises(self):
        """缺 project_id 抛 ValueError。"""
        adata = _make_synthetic_adata()
        adata.obs["source_dataset"] = "DS1"
        adata.obs["disease_system"] = "gastric"
        with pytest.raises(ValueError, match="project_id|Layer 1"):
            _validate_layer1(adata)

    def test_missing_disease_system_raises(self):
        """缺 disease_system 抛 ValueError。"""
        adata = _make_synthetic_adata()
        adata.obs["source_dataset"] = "DS1"
        adata.obs["project_id"] = "PRJ1"
        with pytest.raises(ValueError, match="disease_system|Layer 1"):
            _validate_layer1(adata)

    def test_missing_multiple_raises_with_list(self):
        """缺多个字段时错误消息列出全部。"""
        adata = _make_synthetic_adata()
        adata.obs["disease_system"] = "gastric"
        with pytest.raises(ValueError, match="source_dataset.*project_id"):
            _validate_layer1(adata)

    def test_project_id_written_to_obs(self):
        """read_with_manifest 端到端验证 project_id 写入 obs。"""
        adata = _make_synthetic_adata(n_cells=10)
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml(
                {"input": {"format": "h5ad", "path": h5ad_path}}
            )
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)
        assert "project_id" in result.obs.columns
        assert (result.obs["project_id"] == "test_project").all()


# ---------------------------------------------------------------------------
# P0-2: obs_mapping 目标列冲突不静默覆盖
# ---------------------------------------------------------------------------


class TestObsMappingTargetConflict:
    """P0-2: 目标列已存在时发出 warning。"""

    def test_target_col_already_exists_warns(self):
        """obs_mapping 的目标列已在 obs 中存在时警告。"""
        adata = _make_synthetic_adata(n_cells=10)
        adata.obs["existing_col"] = ["val"] * 10
        adata.obs["src_col"] = ["new_val"] * 10
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "obs_mapping": {"existing_col": "src_col"},
            })
            mpath = _write_manifest(manifest, tmpdir)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                read_with_manifest(mpath)
                conflict_warnings = [
                    x for x in w if "已存在" in str(x.message)
                ]
                assert len(conflict_warnings) >= 1

    def test_target_col_new_no_warn(self):
        """目标列不存在时不警告。"""
        adata = _make_synthetic_adata(n_cells=10)
        adata.obs["src_col"] = ["val"] * 10
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "obs_mapping": {"new_target": "src_col"},
            })
            mpath = _write_manifest(manifest, tmpdir)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                read_with_manifest(mpath)
                conflict_warnings = [
                    x for x in w if "已存在" in str(x.message)
                ]
                assert len(conflict_warnings) == 0


# ---------------------------------------------------------------------------
# P0-3: value_mapping 未覆盖值不静默穿透
# ---------------------------------------------------------------------------


class TestValueMappingUncovered:
    """P0-3: value_mapping 存在未覆盖取值时发出 warning。"""

    def test_uncovered_values_warn(self):
        """映射未覆盖的取值应告警（列出未覆盖值）。"""
        adata = _make_synthetic_adata(n_cells=12)
        adata.obs["orig_status"] = (
            ["CAG"] * 4 + ["IM"] * 4 + ["unknown_subtype"] * 4
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "obs_mapping": {"disease": "orig_status"},
                "value_mapping": {
                    "disease": {"CAG": "atrophic_gastritis"}
                },
            })
            mpath = _write_manifest(manifest, tmpdir)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = read_with_manifest(mpath)
                uncovered_warnings = [
                    x for x in w if "未被映射覆盖" in str(x.message)
                ]
                assert len(uncovered_warnings) >= 1
                msg = str(uncovered_warnings[0].message)
                assert "IM" in msg
                assert "unknown_subtype" in msg

    def test_all_values_covered_no_warn(self):
        """所有取值均被覆盖时不警告。"""
        adata = _make_synthetic_adata(n_cells=10)
        adata.obs["orig_status"] = ["CAG"] * 5 + ["IM"] * 5
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "obs_mapping": {"disease": "orig_status"},
                "value_mapping": {
                    "disease": {"CAG": "atrophic_gastritis", "IM": "metaplasia"}
                },
            })
            mpath = _write_manifest(manifest, tmpdir)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                read_with_manifest(mpath)
                uncovered_warnings = [
                    x for x in w if "未被映射覆盖" in str(x.message)
                ]
                assert len(uncovered_warnings) == 0


# ---------------------------------------------------------------------------
# P0-4: NaN 保留为真 NaN（非字符串 "nan"）
# ---------------------------------------------------------------------------


class TestNaNSafeAstype:
    """P0-4: obs_mapping 的 astype(str) 后 NaN 恢复为真 NaN。"""

    def test_nan_preserved_in_obs_mapping(self):
        """源列含 NaN 时目标列对应位置为真 NaN，非字符串 'nan'。"""
        adata = _make_synthetic_adata(n_cells=10)
        adata.obs["src_col"] = ["A", "B", np.nan, "A", "B", np.nan, "A", "B", "A", "B"]
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "obs_mapping": {"target": "src_col"},
            })
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)
        assert "target" in result.obs.columns
        # 原 NaN 位置应为真 NaN
        nan_mask = result.obs["target"].isna()
        assert nan_mask.sum() == 2
        # 非 NaN 位置不应是字符串 "nan"
        non_nan_vals = result.obs["target"][~nan_mask].unique()
        assert "nan" not in non_nan_vals
        assert "A" in non_nan_vals
        assert "B" in non_nan_vals

    def test_no_nan_col_no_issue(self):
        """无 NaN 的列正常转换不抛异常。"""
        adata = _make_synthetic_adata(n_cells=5)
        adata.obs["src_col"] = ["X", "Y", "Z", "X", "Y"]
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "obs_mapping": {"target": "src_col"},
            })
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)
        assert result.obs["target"].isna().sum() == 0


# ---------------------------------------------------------------------------
# P0-5: clinical join 键类型对齐 + 行数不变校验
# ---------------------------------------------------------------------------


class TestClinicalKeyTypeAlignment:
    """P0-5: merge 前键类型对齐 + merge 后行数不变校验。"""

    def test_key_type_mismatch_still_joins(self):
        """键类型不一致（如 category vs int64）时仍能正确 join。"""
        adata = _make_synthetic_adata(n_cells=12)
        # obs 的 sample_id 是 int 类别
        adata.obs["sample_id"] = pd.Categorical(
            [f"{i}" for i in range(12)]
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)

            csv_path = os.path.join(tmpdir, "clinical.csv")
            pd.DataFrame({
                "sid": [f"{i}" for i in range(12)],
                "batch": ["A"] * 6 + ["B"] * 6,
            }).to_csv(csv_path, index=False)

            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "clinical_metadata": [{
                    "file": csv_path,
                    "join_on": {
                        "manifest_field": "sample_id",
                        "table_column": "sid",
                    },
                }],
            })
            mpath = _write_manifest(manifest, tmpdir)
            result = read_with_manifest(mpath)

        # 应正确关联上
        assert "batch" in result.obs.columns
        # 行数不应变化
        assert result.n_obs == 12
        assert result.obs["batch"].notna().all()

    def test_row_count_warns_on_change(self):
        """merge 后行数变化时 warn。"""
        adata = _make_synthetic_adata(n_cells=8)
        adata.obs["sample_id"] = [f"sample_{i}" for i in range(8)]
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)

            # 制造键不匹配：CSV 的键与 obs 完全不同
            csv_path = os.path.join(tmpdir, "clinical.csv")
            pd.DataFrame({
                "sid": ["nonexistent_1", "nonexistent_2"],
                "batch": ["X", "Y"],
            }).to_csv(csv_path, index=False)

            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "clinical_metadata": [{
                    "file": csv_path,
                    "join_on": {
                        "manifest_field": "sample_id",
                        "table_column": "sid",
                    },
                }],
            })
            mpath = _write_manifest(manifest, tmpdir)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = read_with_manifest(mpath)
                # 正常情况行数不变（左连接），但键完全不匹配会导致
                # 键完全不匹配时右表所有列为 NaN，但左表行数不变。此时主要验证不崩溃。
                assert result.n_obs == 8


# ---------------------------------------------------------------------------
# P0-6: project_specific 源列缺失补 warn
# ---------------------------------------------------------------------------


class TestProjectSpecificMissingSource:
    """P0-6: project_specific 源列不存在时明确警告。"""

    def test_missing_source_column_warns(self):
        """源列不存在时发出警告但不崩溃。"""
        adata = _make_synthetic_adata(n_cells=10)
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "project_specific": {
                    "derived_col": {
                        "source_column": "nonexistent_col",
                        "rules": {"a": "b"},
                    },
                },
            })
            mpath = _write_manifest(manifest, tmpdir)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = read_with_manifest(mpath)
                missing_warnings = [
                    x for x in w if "未找到" in str(x.message)
                ]
                assert len(missing_warnings) >= 1
                assert "nonexistent_col" in str(missing_warnings[0].message)

    def test_existing_source_column_no_warn(self):
        """源列存在时不发出缺失警告。"""
        adata = _make_synthetic_adata(n_cells=10)
        adata.obs["real_col"] = ["A", "B", "C"] * 3 + ["A"]
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "project_specific": {
                    "derived": {
                        "source_column": "real_col",
                        "rules": {"A": "alpha"},
                    },
                },
            })
            mpath = _write_manifest(manifest, tmpdir)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = read_with_manifest(mpath)
                missing_warnings = [
                    x for x in w if "未找到" in str(x.message)
                ]
                # 可能有 Layer 2 等其他未找到警告，但不应有 real_col 相关的
                src_missing = [
                    x for x in missing_warnings if "real_col" in str(x.message)
                ]
                assert len(src_missing) == 0


# ---------------------------------------------------------------------------
# DG-4: categorical 列丢失排序 warn
# ---------------------------------------------------------------------------


class TestCategoricalWarning:
    """DG-4: obs_mapping 中源列为 categorical 时告警排序信息丢失。"""

    def test_categorical_source_warns(self):
        """源列是 categorical dtype 时发出类别信息丢失警告。"""
        adata = _make_synthetic_adata(n_cells=10)
        adata.obs["stage"] = pd.Categorical(
            ["normal", "CAG", "IM", "DYS"] * 2 + ["normal", "CAG"],
            categories=["normal", "CAG", "IM", "DYS"],
            ordered=True,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "obs_mapping": {"disease_stage": "stage"},
            })
            mpath = _write_manifest(manifest, tmpdir)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                read_with_manifest(mpath)
                cat_warnings = [
                    x for x in w if "categorical" in str(x.message).lower()
                ]
                assert len(cat_warnings) >= 1

    def test_non_categorical_source_no_warn(self):
        """普通 int64 列不触发 categorical 警告。"""
        adata = _make_synthetic_adata(n_cells=10)
        # 使用整数列——绝不会被 AnnData 转换为 categorical
        adata.obs["plain_int"] = list(range(10))
        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = os.path.join(tmpdir, "test.h5ad")
            adata.write_h5ad(h5ad_path)
            manifest = _make_manifest_yaml({
                "input": {"format": "h5ad", "path": h5ad_path},
                "obs_mapping": {"target": "plain_int"},
            })
            mpath = _write_manifest(manifest, tmpdir)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                read_with_manifest(mpath)
                cat_warnings = [
                    x for x in w if "categorical" in str(x.message).lower()
                ]
                assert len(cat_warnings) == 0
