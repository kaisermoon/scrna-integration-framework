#!/usr/bin/env python3
"""
make_test_subset.py

Sample representative test fixture subsets from full GCPL datasets (~772k cells).

Produces two fixture sets into `data/_subset/` (gitignored):
  Fixture A: Per-source subsets preserving original data shapes for 01-02 testing.
  Fixture B: Downstream subset from 02_qc_filtered_data.h5ad for 03+ testing.

All random sampling uses a fixed seed for reproducibility.
Each source is sampled independently; one source failure does not crash the whole script.

Usage:
  conda activate scrna-integration
  python scripts/make_test_subset.py
"""

import gc
import gzip
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════════════════════
# PARAMS
# ═══════════════════════════════════════════════════════════════════════════════

SEED = 42
KIM_TARGET = 1500
NANCANG_TARGET = 2000
NANCANG_PER_SAMPLE = 250       # ~cells per selected sample
NANCANG_SAMPLES_PER_GROUP = 2  # per GC/GS/IM
NOWICKI_TARGET = 2500
NOWICKI_MIN_PER_TYPE = 50
TSUBOSAKA_TARGET = 2000
YUE_TARGET = 1000
FIXTURE_B_TARGET = 6000

# Source paths (READ-ONLY, must not write). Based on home directory to
# avoid hard-coding a real username path in the public repository.
GCPL_ROOT = Path.home() / "Works" / "GCPL_scRNA"
KIM_DIR = GCPL_ROOT / "data" / "Kim_2025_Series GSE295401"
NANCANG_DIR = GCPL_ROOT / "data" / "Nancang_2025_GSE249874"
NOWICKI_PATH = str(GCPL_ROOT / "data" / "Nowicki_Osuch_et_al_2023" / "Nowicki_Osuch_et_al_2023.h5ad")
TSUBOSAKA_PATH = str(GCPL_ROOT / "data" / "Tsubosaka_2023_HCA" / "data_9_9_annotated_seurat_all_ut.rds")
YUE_DIR = GCPL_ROOT / "data" / "Yue_SSK_2025_GSE210991"
FIXTURE_B_SOURCE = str(GCPL_ROOT / "results" / "data_objects" / "02_qc_filtered_data.h5ad")

# Output directory (gitignored, within worktree)
OUTPUT_DIR = Path("data/_subset")

# ═══════════════════════════════════════════════════════════════════════════════
# Utilities
# ═══════════════════════════════════════════════════════════════════════════════

rng = np.random.default_rng(SEED)


def _parse_kim_condition(filename: str) -> str:
    """Extract condition from Kim filename.

    Examples:
        GSM8947399_na_GF_1 → na
        GSM8947401_IncomGF_3 → Incom
        GSM8947404_ComBMP_6 → Com
        GSM8947405_SI_GF_9 → SI
        GSM8947407_CN_7 → CN
    """
    parts = Path(filename).stem.split("_")
    if len(parts) < 2:
        return "unknown"
    cond_raw = parts[1]
    # strip treatment suffixes
    for prefix in ("na", "Incom", "Com", "SI", "CN"):
        if cond_raw.startswith(prefix):
            return prefix
    return cond_raw


def _parse_yue_condition(filename: str) -> str:
    """Extract disease condition from Yue filename.

    Examples:
        GSM6443816_GX001-IM3O_count.txt.gz → IM
        GSM6443819_GX023-BO_count.txt.gz → BO
        GSM6443821_GX035-AO_count.txt.gz → AO
    """
    stem = Path(filename).stem.replace("_count", "")
    for part in stem.split("_"):
        if "-" in part:
            sample_id = part.split("-")[1] if "-" in part else part
            for cond in ("IM", "BO", "AO", "AgO", "AL1c", "Nile1O"):
                if sample_id.startswith(cond):
                    # normalize IM3O/IMO → IM, everything else keep as-is
                    if cond in ("IM", "BO", "AO", "AgO"):
                        return cond
                    return sample_id[:3] if len(sample_id) >= 3 else sample_id
    return "unknown"


def _print_distribution(label: str, counts: pd.Series) -> None:
    """Print a compact distribution table."""
    total = counts.sum()
    for k, v in sorted(counts.items()):
        pct = 100 * v / total if total > 0 else 0
        print(f"    {k:30s} {v:6d}  ({pct:5.1f}%)")
    print(f"    {'TOTAL':30s} {total:6d}")


def _save_mtx_10x(adata, out_dir: Path, prefix: str = "") -> None:
    """Write AnnData as 10x-format mtx (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz).

    If prefix is non-empty, writes to out_dir/{prefix}_matrix.mtx.gz etc.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    if prefix:
        mtx_path = out_dir / f"{prefix}_matrix.mtx.gz"
        bc_path = out_dir / f"{prefix}_barcodes.tsv.gz"
        ft_path = out_dir / f"{prefix}_features.tsv.gz"
    else:
        mtx_path = out_dir / "matrix.mtx.gz"
        bc_path = out_dir / "barcodes.tsv.gz"
        ft_path = out_dir / "features.tsv.gz"

    # Write sparse matrix
    with gzip.open(mtx_path, "wb") as f:
        scipy.io.mmwrite(f, adata.X.T)  # 10x convention: genes × cells

    # Write barcodes
    with gzip.open(bc_path, "wt") as f:
        f.write("\n".join(adata.obs_names.tolist()) + "\n")

    # Write features (gene_id + gene_name + feature_type)
    with gzip.open(ft_path, "wt") as f:
        for g in adata.var_names:
            f.write(f"{g}\t{g}\tGene Expression\n")


def _index_sample(rng: np.random.Generator, n_total: int, n_target: int) -> np.ndarray:
    """Return random indices of size n_target (or all if n_total < n_target)."""
    n = min(n_total, n_target)
    return rng.choice(n_total, size=n, replace=False)


# ═══════════════════════════════════════════════════════════════════════════════
# Kim (10x .h5 files, ~1500 cells, stratified by condition)
# ═══════════════════════════════════════════════════════════════════════════════

def sample_kim() -> bool:
    print("=" * 60)
    print("KIM — 10x .h5 files (na/Incom/Com/CN/SI)")
    print("=" * 60)

    h5_files = sorted(KIM_DIR.glob("*.h5"))
    if not h5_files:
        print("  FAIL: no .h5 files found in", KIM_DIR)
        return False

    adatas = []
    for fp in h5_files:
        try:
            adata = sc.read_10x_h5(fp)
            adata.var_names_make_unique()
            cond = _parse_kim_condition(fp.name)
            adata.obs["condition"] = cond
            adata.obs["source_file"] = fp.name
            adatas.append(adata)
            print(f"  Loaded {fp.name}: {adata.n_obs} cells, condition={cond}")
        except Exception as e:
            print(f"  WARN: failed to load {fp.name}: {e}")
    if not adatas:
        print("  FAIL: could not load any Kim file")
        return False

    merged = adatas[0].concatenate(adatas[1:], batch_key="_batch", join="outer")
    print(f"  Merged: {merged.n_obs} cells, {merged.n_vars} genes")

    # Stratified sample by condition
    conditions = merged.obs["condition"].value_counts()
    print("  Condition distribution (full):")
    _print_distribution("", conditions)

    frac = KIM_TARGET / merged.n_obs
    sampled = []
    for cond in conditions.index:
        mask = merged.obs["condition"] == cond
        idx = np.where(mask)[0]
        n_take = max(1, int(len(idx) * frac))
        take = rng.choice(idx, size=min(n_take, len(idx)), replace=False)
        sampled.append(take)
    indices = np.sort(np.concatenate(sampled))
    adata_out = merged[indices].copy()

    print(f"  Sampled: {adata_out.n_obs} cells")
    print("  Condition distribution (sampled):")
    _print_distribution("", adata_out.obs["condition"].value_counts())

    out_dir = OUTPUT_DIR / "kim"
    out_dir.mkdir(parents=True, exist_ok=True)
    adata_out.write_h5ad(out_dir / "kim_subset.h5ad", compression="lzf")
    print(f"  Wrote {out_dir / 'kim_subset.h5ad'}")
    del adata_out, merged, adatas
    gc.collect()
    return True


# ═══════════════════════════════════════════════════════════════════════════════
# Nancang (10x mtx, ~2000 cells, filtered+raw, stratified GC/GS/IM)
# ═══════════════════════════════════════════════════════════════════════════════

def sample_nancang() -> bool:
    print("\n" + "=" * 60)
    print("NANCANG — 10x mtx (GC/GS/IM, filtered + raw)")
    print("=" * 60)

    # Disease group → GSM accessions (from metadata)
    groups = {
        "GC": ["GSM7966226", "GSM7966227", "GSM7966228", "GSM7966229", "GSM7966230", "GSM7966231"],
        "GS": ["GSM7966232", "GSM7966233", "GSM7966234", "GSM7966235", "GSM7966236", "GSM7966237"],
        "IM": ["GSM7966238", "GSM7966239", "GSM7966240", "GSM7966241", "GSM7966242", "GSM7966243"],
    }

    total_sampled = 0
    for group_name, accessions in groups.items():
        selected = accessions[:NANCANG_SAMPLES_PER_GROUP]
        for acc in selected:
            src_filtered = NANCANG_DIR / acc / "filtered_feature_bc_matrix"
            src_raw = NANCANG_DIR / acc / "raw_feature_bc_matrix"
            if not (src_filtered / "matrix.mtx.gz").exists():
                print(f"  WARN: {acc} filtered mtx not found, skipping")
                continue

            try:
                adata = sc.read_10x_mtx(src_filtered, var_names="gene_symbols")
                adata.var_names_make_unique()
                print(f"  Loaded {acc} ({group_name}): {adata.n_obs} cells (filtered)")

                # Sample cells
                n_sample = min(NANCANG_PER_SAMPLE, adata.n_obs)
                idx = rng.choice(adata.n_obs, size=n_sample, replace=False)
                barcodes = adata.obs_names[idx]

                # Write filtered subset as 10x mtx
                adata_filt = adata[idx].copy()
                out_filtered = OUTPUT_DIR / "nancang" / acc / "filtered_feature_bc_matrix"
                _save_mtx_10x(adata_filt, out_filtered)

                # Sample matching cells from raw matrix
                raw_exists = (src_raw / "matrix.mtx.gz").exists()
                if raw_exists:
                    adata_raw = sc.read_10x_mtx(src_raw, var_names="gene_symbols")
                    adata_raw.var_names_make_unique()
                    # Match barcodes (may have suffix like '-1')
                    raw_barcodes = set(adata_raw.obs_names)
                    matched = [b for b in barcodes if b in raw_barcodes]
                    if matched:
                        adata_raw_sub = adata_raw[matched].copy()
                        out_raw = OUTPUT_DIR / "nancang" / acc / "raw_feature_bc_matrix"
                        _save_mtx_10x(adata_raw_sub, out_raw)
                        print(f"    Raw: {len(matched)} matching barcodes written")
                    else:
                        print("    Raw: no matching barcodes found (skipping raw)")
                        del adata_raw
                else:
                    print(f"    Raw: not available for {acc}")

                total_sampled += n_sample
                del adata, adata_filt
                gc.collect()

            except Exception as e:
                print(f"  FAIL: {acc}: {e}")
                continue

    print(f"  Total sampled across all Nancang samples: {total_sampled}")
    return total_sampled > 0


# ═══════════════════════════════════════════════════════════════════════════════
# Nowicki (h5ad, ~2500 cells, stratified by Celltypes_global × Patient_status)
# ═══════════════════════════════════════════════════════════════════════════════

def sample_nowicki() -> bool:
    print("\n" + "=" * 60)
    print("NOWICKI — h5ad (Celltypes_global × Patient_status stratified)")
    print("=" * 60)

    if not os.path.exists(NOWICKI_PATH):
        print(f"  FAIL: file not found: {NOWICKI_PATH}")
        return False

    try:
        adata = sc.read_h5ad(NOWICKI_PATH)
        print(f"  Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    except Exception as e:
        print(f"  FAIL: cannot read {NOWICKI_PATH}: {e}")
        return False

    if "Celltypes_global" not in adata.obs.columns:
        print("  FAIL: 'Celltypes_global' column not found in obs")
        return False
    if "Patient_status" not in adata.obs.columns:
        print("  FAIL: 'Patient_status' column not found in obs")
        return False

    ct_col = "Celltypes_global"
    ps_col = "Patient_status"
    cross = adata.obs.groupby([ct_col, ps_col]).size()
    print(f"  Cross-tab groups: {len(cross)}")
    n_total = adata.n_obs

    # Stratified sampling: proportional to group size, floor at MIN_PER_TYPE per CT
    indices = []
    for (ct, ps), grp_n in cross.items():
        # Determine target for this (ct, ps) group
        prop = grp_n / n_total
        n_from_prop = int(NOWICKI_TARGET * prop)
        n_take = max(NOWICKI_MIN_PER_TYPE, n_from_prop) if ct in set(
            cross.index.get_level_values(0)[cross >= NOWICKI_MIN_PER_TYPE]
        ) else n_from_prop
        n_take = max(1, n_take)
        mask = (adata.obs[ct_col] == ct) & (adata.obs[ps_col] == ps)
        grp_idx = np.where(mask)[0]
        n_take = min(n_take, len(grp_idx))
        take = rng.choice(grp_idx, size=n_take, replace=False)
        indices.append(take)

    all_idx = np.sort(np.concatenate(indices))
    # If we overshot, randomly trim
    if len(all_idx) > NOWICKI_TARGET:
        all_idx = rng.choice(all_idx, size=NOWICKI_TARGET, replace=False)
        all_idx.sort()

    adata_out = adata[all_idx].copy()
    print(f"  Sampled: {adata_out.n_obs} cells (target: {NOWICKI_TARGET})")
    print("  Celltypes_global distribution:")
    _print_distribution("", adata_out.obs[ct_col].value_counts())
    print("  Patient_status distribution:")
    _print_distribution("", adata_out.obs[ps_col].value_counts())

    out_dir = OUTPUT_DIR / "nowicki"
    out_dir.mkdir(parents=True, exist_ok=True)
    adata_out.write_h5ad(out_dir / "nowicki_subset.h5ad", compression="lzf")
    print(f"  Wrote {out_dir / 'nowicki_subset.h5ad'}")

    del adata, adata_out
    gc.collect()
    return True


# ═══════════════════════════════════════════════════════════════════════════════
# Tsubosaka (Seurat RDS) — TODO
# ═══════════════════════════════════════════════════════════════════════════════

def sample_tsubosaka() -> bool:
    print("\n" + "=" * 60)
    print("TSUBOSAKA — Seurat RDS (major_clusters × subtype stratified)")
    print("=" * 60)

    # rds2py (pure Python) cannot parse Seurat S4 objects with deeply nested
    # pairlists (observed error: "failed to parse pairlist element 'commands' /
    # cannot read unknown SEXP type 3"). The RDS stores command history as
    # serialized closures which rds2py does not support.
    #
    # Fallback: R conda environment + subprocess Rscript.
    # Implementation sketch for when the R env is set up:
    #
    #   library(Seurat)
    #   so <- readRDS("data_9_9_annotated_seurat_all_ut.rds")
    #   # stratified sample by major_clusters × subtype
    #   target <- 2000
    #   set.seed(42)
    #   cells <- colnames(so)
    #   meta <- so@meta.data[, c("major_clusters", "subtype", "subcluster")]
    #   sampled <- ...
    #   so_sub <- so[, sampled]
    #   # Convert to h5ad via SeuratDisk or anndata write
    #   SaveH5Seurat(so_sub, "tsubosaka_subset.h5seurat")
    #   Convert("tsubosaka_subset.h5seurat", dest="h5ad")
    #
    # For now, skip with a clear message.

    print("  SKIPPED: rds2py cannot parse Seurat S4 with embedded command")
    print("           history. Requires R conda environment + subprocess Rscript.")
    print("           See comments in sample_tsubosaka() for the R fallback sketch.")
    print("           TARGET: ~2000 cells stratified by major_clusters × subtype.")

    # Create an empty marker so the output dir exists
    out_dir = OUTPUT_DIR / "tsubosaka"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "_TODO_R_setup_needed.txt", "w") as f:
        f.write(
            "Tsubosaka RDS sampling requires R environment.\n"
            "See scripts/make_test_subset.py sample_tsubosaka() for the Rscript fallback.\n"
        )
    return False  # not a crash — explicitly skipped


# ═══════════════════════════════════════════════════════════════════════════════
# Yue (txt.gz counts, ~1000 cells, covering IM/BO/AO)
# ═══════════════════════════════════════════════════════════════════════════════

def sample_yue() -> bool:
    print("\n" + "=" * 60)
    print("YUE — txt.gz counts (IM/BO/AO)")
    print("=" * 60)

    # Only *_count.txt.gz files, skip *_LogNormalized.txt.gz
    count_files = sorted(YUE_DIR.glob("*_count.txt.gz"))
    if not count_files:
        print("  FAIL: no count.txt.gz files found in", YUE_DIR)
        return False

    # Group by condition
    file_groups: dict[str, list[Path]] = {}
    for fp in count_files:
        cond = _parse_yue_condition(fp.name)
        file_groups.setdefault(cond, []).append(fp)
    print(f"  Found conditions: {list(file_groups.keys())}")
    for cond, files in sorted(file_groups.items()):
        print(f"    {cond}: {len(files)} files")

    # Target distribution: proportional to file count per condition
    total_files = len(count_files)
    condition_counts: dict[str, int] = {}

    for cond, files in sorted(file_groups.items()):
        cond_target = max(50, int(YUE_TARGET * len(files) / total_files))
        cond_sampled = 0

        for fp in files:
            if cond_sampled >= cond_target:
                break
            try:
                # Read count matrix (genes × cells in typical txt.gz export)
                df = pd.read_csv(fp, sep="\t", index_col=0, compression="gzip")
                df = df.T  # cells × genes
                if df.shape[0] == 0:
                    print(f"    WARN: {fp.name} empty, skipping")
                    continue

                n_take = min(cond_target - cond_sampled, df.shape[0],
                             max(1, cond_target // len(files)))
                idx = rng.choice(df.shape[0], size=n_take, replace=False)
                sampled_df = df.iloc[idx]

                # Write sampled subset back as txt.gz (genes × cells)
                out_dir = OUTPUT_DIR / "yue"
                out_dir.mkdir(parents=True, exist_ok=True)
                out_path = out_dir / fp.name
                sampled_df.T.to_csv(out_path, sep="\t", compression="gzip")
                cond_sampled += n_take
                condition_counts[cond] = condition_counts.get(cond, 0) + n_take
                print(f"    {fp.name}: {n_take} cells → {out_path.name}")

            except Exception as e:
                print(f"    FAIL: {fp.name}: {e}")
                continue

    print(f"  Total sampled: {sum(condition_counts.values())}")
    print("  Condition distribution:")
    for cond, n in sorted(condition_counts.items()):
        print(f"    {cond:10s} {n:6d}")

    gc.collect()
    return sum(condition_counts.values()) > 0


# ═══════════════════════════════════════════════════════════════════════════════
# Fixture B — downstream subset from 02_qc_filtered_data.h5ad (6.9G)
# ═══════════════════════════════════════════════════════════════════════════════

def sample_fixture_b() -> bool:
    print("\n" + "=" * 60)
    print("FIXTURE B — from 02_qc_filtered_data.h5ad (6.9G, backed='r')")
    print("=" * 60)

    if not os.path.exists(FIXTURE_B_SOURCE):
        print(f"  FAIL: file not found: {FIXTURE_B_SOURCE}")
        return False

    try:
        # Use backed='r' to avoid loading 6.9G into memory
        adata = sc.read_h5ad(FIXTURE_B_SOURCE, backed="r")
        n_total = adata.n_obs
        print(f"  Loaded (backed): {n_total} cells, {adata.n_vars} genes")
        print(f"  obs columns: {list(adata.obs.columns)[:10]}...")
    except Exception as e:
        print(f"  FAIL: cannot open with backed='r': {e}")
        return False

    # Determine stratification column (disease-related)
    # Check for common disease column names
    obs_cols = list(adata.obs.columns)
    disease_col = None
    for candidate in ["disease", "group", "condition", "Group", "disease_status",
                      "patient_status", "Patient_status", "sample_group"]:
        if candidate in obs_cols:
            disease_col = candidate
            break

    if disease_col is not None:
        print(f"  Stratifying by: {disease_col}")
        groups = adata.obs[disease_col].value_counts()
        print("  Group distribution (full):")
        _print_distribution("", groups)

        # Proportional stratified sample
        sampled_idx = []
        frac = FIXTURE_B_TARGET / n_total
        for grp in groups.index:
            mask = adata.obs[disease_col] == grp
            grp_indices = np.where(mask)[0]
            n_take = max(1, int(len(grp_indices) * frac))
            n_take = min(n_take, len(grp_indices))
            take = rng.choice(grp_indices, size=n_take, replace=False)
            sampled_idx.append(take)
        indices = np.sort(np.concatenate(sampled_idx))
    else:
        print("  No disease column found, random sampling")
        indices = rng.choice(n_total, size=min(FIXTURE_B_TARGET, n_total), replace=False)
        indices.sort()

    if len(indices) > FIXTURE_B_TARGET + 500:
        indices = rng.choice(indices, size=FIXTURE_B_TARGET, replace=False)
        indices.sort()

    print(f"  Loading {len(indices)} sampled cells into memory...")
    adata_out = adata[indices].to_memory()
    print(f"  Sampled: {adata_out.n_obs} cells")

    if disease_col is not None:
        print(f"  {disease_col} distribution:")
        _print_distribution("", adata_out.obs[disease_col].value_counts())

    out_path = OUTPUT_DIR / "fixture_B_qcd.h5ad"
    adata_out.write_h5ad(out_path, compression="lzf")
    print(f"  Wrote {out_path} ({out_path.stat().st_size / 1e6:.1f} MB)")

    del adata, adata_out
    gc.collect()
    return True


# ═══════════════════════════════════════════════════════════════════════════════
# 10x_h5 夹具：从 kim_subset.h5ad 按 source_file 拆回逐样本 CellRanger v3 h5
# ═══════════════════════════════════════════════════════════════════════════════

def _save_10x_h5(adata_out, out_path: Path) -> None:
    """将 AnnData 写出为 CellRanger v3 格式的 10x h5 文件。

    布局：
      matrix/data, indices, indptr, shape
      matrix/features/id, name, feature_type, genome
      matrix/barcodes
    根节点 attrs：chemistry_description, filetype, version
    """
    import h5py  # 仅在 scrna-integration 环境可用（3.14.0）

    out_path.parent.mkdir(parents=True, exist_ok=True)
    mat = adata_out.X
    if not sp.issparse(mat):
        mat = sp.csr_matrix(mat)
    # 10x h5 存为 genes × cells（与 mtx 格式一致），CSC column-major
    mat_10x = mat.T.tocsc()  # cells×genes → genes×cells, CSC

    with h5py.File(out_path, "w") as f:
        # 根节点 attrs
        f.attrs["chemistry_description"] = "Single Cell 3' v3"
        f.attrs["filetype"] = "h5"
        f.attrs["library_ids"] = ["unknown"]
        f.attrs["original_gem_groups"] = ["unknown"]
        f.attrs["version"] = 3

        # matrix group
        mg = f.create_group("matrix")
        mg.create_dataset("data", data=mat_10x.data, dtype="int32")
        mg.create_dataset("indices", data=mat_10x.indices, dtype="int32")
        mg.create_dataset("indptr", data=mat_10x.indptr, dtype="int64")
        mg.create_dataset("shape", data=np.array([mat_10x.shape[0], mat_10x.shape[1]], dtype="int32"))

        # features group
        fg = mg.create_group("features")
        gene_ids = adata_out.var["gene_ids"].values.astype("S")
        gene_names = adata_out.var_names.values.astype("S")
        ft = adata_out.var["feature_types"].values.astype("S")
        fg.create_dataset("id", data=gene_ids)
        fg.create_dataset("name", data=gene_names)
        fg.create_dataset("feature_type", data=ft)
        fg.create_dataset(
            "genome", data=np.array(["GRCh38"] * adata_out.n_vars, dtype="S")
        )
        fg.create_dataset("_all_tag_keys", data=np.array([b"genome"], dtype="S"))

        # barcodes
        barcodes = adata_out.obs_names.values.astype("S")
        mg.create_dataset("barcodes", data=barcodes)


def sample_10x_h5_fixtures() -> bool:
    """从 data/_subset/kim/kim_subset.h5ad 按 source_file 拆回逐样本 10x h5 文件。

    产出：
      data/_subset/tenx_h5/<sample_id>_filtered_feature_bc_matrix.h5（逐样本 h5）
      data/_subset/tenx_h5/manifest.yaml（附带 gsm_to_sample 占位映射）

    h5py 由 _save_10x_h5() 内部按需导入（仅 scrna-integration 环境可用），此处不重复导入。
    """
    print()
    print("=" * 60)
    print("10x_h5 夹具 — 从 kim_subset.h5ad 按 source_file 拆分")
    print("=" * 60)

    source_h5ad = OUTPUT_DIR / "kim" / "kim_subset.h5ad"
    if not source_h5ad.exists():
        print(f"  SKIP: source not found: {source_h5ad}")
        return False

    adata = sc.read_h5ad(source_h5ad)
    print(f"  源: {adata.n_obs} 细胞 x {adata.n_vars} 基因, {adata.obs['source_file'].nunique()} 个样本")

    out_dir = OUTPUT_DIR / "tenx_h5"
    if out_dir.exists():
        import shutil
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_map = {}  # sample_id → 占位 accession
    for _sf, group in adata.obs.groupby("source_file"):
        # 样本 ID 与 accession 全部用占位命名。
        # 为什么不沿用源文件名里的真实 accession：本仓为 public，夹具 manifest
        # 是模板的参考样例，可能被拷进进 git 的 per-dataset manifest。
        # 真实 accession 一旦入库不可逆，因此在生成阶段就切断，不留侥幸。
        idx = len(sample_map) + 1
        sample_id = f"sample_{idx:02d}"
        placeholder_accession = f"SAMPLE{idx:05d}"
        sample_map[sample_id] = {
            "accession": placeholder_accession,
            "n_cells": len(group),
        }

        # 提取该样本的子集
        ad_sub = adata[group.index].copy()
        # 去掉 source_file 列（含真实文件名，夹具里不需要且不应保留）
        if "source_file" in ad_sub.obs.columns:
            del ad_sub.obs["source_file"]

        # 写出 10x h5（CellRanger v3 布局）
        out_h5 = out_dir / f"{sample_id}_filtered_feature_bc_matrix.h5"
        _save_10x_h5(ad_sub, out_h5)
        print(f"  {sample_id}: {ad_sub.n_obs} 细胞 → {out_h5.name} ({out_h5.stat().st_size / 1e6:.1f} MB)")
        del ad_sub
        gc.collect()

    # 写出 manifest.yaml（含 gsm_to_sample 占位映射）
    import yaml
    gsm_to_sample = {v["accession"]: k for k, v in sample_map.items()}
    manifest = {
        "source_dataset": "TEMPLATE_dataset",
        "project_id": "TEMPLATE_project",
        "disease_system": "TEMPLATE",
        "input": {
            "format": "10x_h5",
            "path": "data/_subset/tenx_h5",
        },
        "gsm_to_sample": gsm_to_sample,  # 占位：GSM 号→样本名映射
        "obs_mapping": {},  # 此格式无内置 obs，下游按需映射
        "preprocessing_done": [],
        "qc_overrides": {},
    }
    manifest_path = out_dir / "manifest.yaml"
    with open(manifest_path, "w") as f:
        yaml.dump(manifest, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
    print(f"  manifest: {manifest_path}")
    print(f"  gsm_to_sample 映射 ({len(gsm_to_sample)} 项)")

    del adata
    gc.collect()
    return True


def _write_format_manifest(out_dir: Path, fmt: str, sample_ids: list[str], **extra) -> Path:
    """写出按格式命名的夹具 manifest.yaml（字段与 tenx_h5 那份保持一致）。

    为什么 source_dataset / project_id / disease_system 全用 TEMPLATE_* 占位：
    本仓为 public，这份 manifest 是模板的参考样例，很可能被拷进下游进 git 的
    per-dataset manifest。真实数据集名、内部课题名、疾病方向一旦入库不可逆，
    因此在生成阶段就切断，不留侥幸。
    """
    import yaml

    manifest = {
        "source_dataset": "TEMPLATE_dataset",
        "project_id": "TEMPLATE_project",
        "disease_system": "TEMPLATE",
        "input": {"format": fmt, "path": f"data/_subset/{out_dir.name}"},
        # 占位 accession→样本名映射：真实项目里 key 换成 GEO 的 Sample_title 对应编号
        "gsm_to_sample": {
            f"SAMPLE{i:05d}": sid for i, sid in enumerate(sample_ids, start=1)
        },
        "obs_mapping": {},
        "preprocessing_done": [],
        "qc_overrides": {},
    }
    manifest.update(extra)
    manifest_path = out_dir / "manifest.yaml"
    with open(manifest_path, "w") as f:
        yaml.dump(manifest, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
    print(f"  manifest: {manifest_path}（{len(sample_ids)} 个样本占位映射）")
    return manifest_path


def _reset_dir(out_dir: Path) -> None:
    """清空并重建输出目录，保证夹具可重复生成、不残留上一次的文件。"""
    if out_dir.exists():
        import shutil

        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)


def sample_10x_mtx_fixtures() -> bool:
    """从 data/_subset/nancang/ 转写出按格式命名的 10x mtx 夹具。

    产出：
      data/_subset/tenx_mtx/sample_NN/filtered_feature_bc_matrix/{matrix.mtx.gz,...}
      data/_subset/tenx_mtx/sample_NN/raw_feature_bc_matrix/（若源目录有）
      data/_subset/tenx_mtx/manifest.yaml

    为什么必须保留 filtered/raw 两层目录结构而不能扁平化：
    01_template_10x_mtx 的子目录发现逻辑优先识别 CellRanger 双层结构
    （sample_NN/filtered_feature_bc_matrix/），并据此推导 SoupX 需要的
    raw droplets 路径（同级的 raw_feature_bc_matrix/）。拍平成单层后
    SoupX 找不到 raw 矩阵，环境 RNA 校正整段会被跳过。
    """
    import shutil

    print()
    print("=" * 60)
    print("10x_mtx 夹具 — 从按数据集命名的 mtx 夹具转写")
    print("=" * 60)

    src_root = OUTPUT_DIR / "nancang"
    if not src_root.is_dir():
        print(f"  SKIP: source not found: {src_root}")
        return False

    # 源目录下每个子目录是一个样本（目录名含真实 accession，转写时丢弃）
    src_samples = sorted(d for d in src_root.iterdir() if d.is_dir())
    if not src_samples:
        print(f"  SKIP: no sample subdirectories under {src_root}")
        return False

    out_dir = OUTPUT_DIR / "tenx_mtx"
    _reset_dir(out_dir)

    sample_ids = []
    for idx, src_sample in enumerate(src_samples, start=1):
        # 目录名一律改占位：源目录名含真实 GEO accession，不能带进夹具
        sample_id = f"sample_{idx:02d}"
        sample_ids.append(sample_id)
        dst_sample = out_dir / sample_id
        dst_sample.mkdir(parents=True, exist_ok=True)

        copied = []
        for layer in ("filtered_feature_bc_matrix", "raw_feature_bc_matrix"):
            src_layer = src_sample / layer
            if src_layer.is_dir():
                shutil.copytree(src_layer, dst_sample / layer)
                copied.append(layer)
        if not copied:
            # 源是扁平结构：包一层 filtered_feature_bc_matrix，对齐模板期望
            dst_layer = dst_sample / "filtered_feature_bc_matrix"
            dst_layer.mkdir(parents=True, exist_ok=True)
            for f in src_sample.iterdir():
                if f.is_file():
                    shutil.copy2(f, dst_layer / f.name)
            copied.append("filtered_feature_bc_matrix(由扁平结构包装)")
        print(f"  {sample_id}: {', '.join(copied)}")

    _write_format_manifest(out_dir, "10x_mtx", sample_ids)
    return True


def sample_h5ad_fixtures() -> bool:
    """从按数据集命名的预处理 h5ad 夹具转写出按格式命名的 h5ad 夹具。

    产出：
      data/_subset/h5ad/sample_subset.h5ad
      data/_subset/h5ad/manifest.yaml（preprocessing_done + qc_overrides 声明作者已 QC）

    为什么必须保留 .raw：
    01_template_h5ad 的核心路径是从 .raw.X 把原始 counts 恢复到 layers["counts"]
    （预处理 h5ad 的 X 通常已是 log-normalized，不能当 counts 用）。丢掉 .raw
    该模板的 expression_contract 建立步骤就没有数据来源。
    """
    print()
    print("=" * 60)
    print("h5ad 夹具 — 从按数据集命名的预处理 h5ad 转写")
    print("=" * 60)

    src = OUTPUT_DIR / "nowicki" / "nowicki_subset.h5ad"
    if not src.exists():
        print(f"  SKIP: source not found: {src}")
        return False

    adata = sc.read_h5ad(src)
    print(f"  源: {adata.n_obs} 细胞 x {adata.n_vars} 基因, .raw={'有' if adata.raw is not None else '无'}")

    out_dir = OUTPUT_DIR / "h5ad"
    _reset_dir(out_dir)

    # 脱敏 obs 中可能含真实标识的列。
    # nowicki 源数据的标识列：Sample(GSM accession)、donor_id(患者ID)、
    # Study(真实研究名)、Sample_Barcode(条码)、Batch(批次号)。
    # 为什么要脱敏：本仓为 public，夹具文件虽不进 git 但可能被分享，
    # 任何可追溯到真实患者/研究的标识必须在生成阶段切断。
    sample_ids = ["sample_01"]  # h5ad 是单文件夹具，manifest 用单样本占位
    for col, placeholder_prefix in [
        ("Sample", "SAMPLE"),
        ("donor_id", "DONOR"),
        ("Sample_Barcode", "BARCODE"),
        ("Batch", "BATCH"),
        ("Study", "TEMPLATE_study"),
    ]:
        if col not in adata.obs.columns:
            continue
        # 先转为 object dtype（保留 NaN），避免 Categorical setitem 报错；
        # astype(str) 会把 NaN 转成字面串 "nan"，所以必须在转换前存好 NaN 掩码
        mask_notna = adata.obs[col].notna()
        adata.obs[col] = adata.obs[col].astype(object)
        # 取所有非 NaN 的非占位唯一值（"NA" 已是占位，无需脱敏）
        vals = sorted(
            str(v) for v in adata.obs.loc[mask_notna, col].unique()
            if str(v) != "NA"
        )
        if not vals:
            print(f"  obs['{col}']: 无需脱敏（全为 NA/nan）")
            continue
        mapping = {v: f"{placeholder_prefix}{i:05d}" for i, v in enumerate(vals, start=1)}
        adata.obs.loc[mask_notna, col] = (
            adata.obs.loc[mask_notna, col].astype(str).replace(mapping)
        )
        print(f"  obs['{col}']: {len(mapping)} 个取值已脱敏")

    out_h5ad = out_dir / "sample_subset.h5ad"
    adata.write_h5ad(out_h5ad, compression="gzip")
    print(f"  写出: {out_h5ad.name} ({out_h5ad.stat().st_size / 1e6:.1f} MB)")

    # 本格式的关键语义：原作者已完成 QC + normalization，重复过滤是科学错误，
    # 故 manifest 用 qc_overrides 逐步声明 skip，模板的 preflight 据此走跳过分支。
    done = ["basic_filter", "doublet_removal", "normalization"]
    _write_format_manifest(
        out_dir, "h5ad", sample_ids,
        preprocessing_done=done,
        qc_overrides={s: {"skip": True, "reason": "原作者已完成"} for s in done},
    )
    # input.path 指向具体文件而非目录（本格式是单文件输入）
    import yaml

    mpath = out_dir / "manifest.yaml"
    with open(mpath) as f:
        m = yaml.safe_load(f)
    m["input"]["path"] = f"data/_subset/h5ad/{out_h5ad.name}"
    m["input"]["gene_id_format"] = "ensembl"  # 预处理 h5ad 常见为 Ensembl id
    with open(mpath, "w") as f:
        yaml.dump(m, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

    del adata
    gc.collect()
    return True


def sample_counts_matrix_fixtures() -> bool:
    """从按数据集命名的 gzip 计数表夹具转写出按格式命名的 counts_matrix 夹具。

    产出：
      data/_subset/counts_matrix/sample_NN_count.txt.gz
      data/_subset/counts_matrix/manifest.yaml

    为什么要重命名文件：源文件名形如 <accession>_<内部样本编号>_count.txt.gz，
    同时含真实 GEO accession 与真实样本编号。本仓为 public，夹具文件名会出现在
    模板的 print 输出与 manifest 里，必须在生成阶段就换成占位名。
    格式（gzip + tab 分隔 + 基因行名）原样保留，不重新编码。
    """

    print()
    print("=" * 60)
    print("counts_matrix 夹具 — 从按数据集命名的 gzip 计数表转写")
    print("=" * 60)

    src_root = OUTPUT_DIR / "yue"
    if not src_root.is_dir():
        print(f"  SKIP: source not found: {src_root}")
        return False

    src_files = sorted(
        f for f in src_root.iterdir()
        if f.is_file() and f.name.endswith((".txt.gz", ".tsv.gz"))
    )
    if not src_files:
        print(f"  SKIP: no .txt.gz/.tsv.gz files under {src_root}")
        return False

    out_dir = OUTPUT_DIR / "counts_matrix"
    _reset_dir(out_dir)

    sample_ids = []
    for idx, src_file in enumerate(src_files, start=1):
        sample_id = f"sample_{idx:02d}"
        sample_ids.append(sample_id)
        dst = out_dir / f"{sample_id}_count.txt.gz"
        # 不能直接用 shutil.copy2：gzip header 的 FNAME 字段存了源文件名
        # （形如 GSM6443816_GX001-IM3O_count.txt），直接复制会把真实 accession
        # 带进夹具文件。解压后重新压缩会写入干净的目标文件名。
        with gzip.open(src_file, "rt") as fin, gzip.open(dst, "wt") as fout:
            fout.write(fin.read())
    print(f"  {len(sample_ids)} 个样本 → sample_01_count.txt.gz … sample_{len(sample_ids):02d}_count.txt.gz")

    _write_format_manifest(out_dir, "tsv_matrix", sample_ids)
    return True


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    print("make_test_subset.py — GCPL fixture sampling")
    print(f"Seed: {SEED}  |  Output: {OUTPUT_DIR.resolve()}")
    print()

    results: dict[str, bool] = {}

    # Fixture A
    for name, fn in [
        ("Kim", sample_kim),
        ("Nancang", sample_nancang),
        ("Nowicki", sample_nowicki),
        ("Tsubosaka", sample_tsubosaka),
        ("Yue", sample_yue),
    ]:
        try:
            results[name] = fn()
        except Exception as e:
            print(f"\n  UNEXPECTED FAILURE in {name}: {e}")
            import traceback
            traceback.print_exc()
            results[name] = False

    # Fixture B
    try:
        results["Fixture_B"] = sample_fixture_b()
    except Exception as e:
        print(f"\n  UNEXPECTED FAILURE in Fixture_B: {e}")
        import traceback
        traceback.print_exc()
        results["Fixture_B"] = False

    # 按格式命名的夹具：供 4 个 01_template_<fmt> 模板直接 Run All。
    # 全部放在最后，因为它们都从上面按数据集命名的夹具转写而来。
    # 为什么要另存一套：模板按输入数据格式命名（数据集属下游项目资产，
    # 格式才是框架该提供的可复用骨架），其 MANIFEST_PATH 默认指向
    # data/_subset/<fmt>/，且文件名与 manifest 全用占位、不含真实 accession。
    for _key, _fn in [
        ("10x_h5", sample_10x_h5_fixtures),
        ("10x_mtx", sample_10x_mtx_fixtures),
        ("h5ad_fmt", sample_h5ad_fixtures),
        ("counts_matrix", sample_counts_matrix_fixtures),
    ]:
        try:
            results[_key] = _fn()
        except Exception as e:
            print(f"\n  UNEXPECTED FAILURE in {_key}: {e}")
            import traceback
            traceback.print_exc()
            results[_key] = False

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    success = 0
    for name, ok in results.items():
        status = "PASS" if ok else ("SKIPPED" if name == "Tsubosaka" else "FAIL")
        print(f"  {name:20s}  {status}")
        if ok or (name == "Tsubosaka" and not ok):
            success += 1
    print(f"\n  {success}/{len(results)} sources completed (Tsubosaka skipped = expected)")

    if all(results.get(k, True) for k in ["Kim", "Nancang", "Nowicki", "Yue", "Fixture_B"]):
        print("  Fixture A (4/5 sources) + Fixture B ready for 01-03+ testing.")
    else:
        print("  Some sources failed unexpectedly — check output above.")


if __name__ == "__main__":
    main()
