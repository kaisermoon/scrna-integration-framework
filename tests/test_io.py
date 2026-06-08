"""Tests for scrna_integration.io — read_with_manifest and retained helpers."""

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
