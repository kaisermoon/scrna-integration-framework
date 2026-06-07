#!/usr/bin/env python3
"""
make_test_subset.py

Sample representative test fixture subsets from full GCPL datasets (~772k cells).

Produces two fixture sets into `data/_subset/` (gitignored):
  Fixture A: Per-source subsets preserving original data shapes for stage1-2 testing.
  Fixture B: Downstream subset from 02_qc_filtered_data.h5ad for stage3+ testing.

All random sampling uses a fixed seed for reproducibility.
Each source is sampled independently; one source failure does not crash the whole script.

Usage:
  conda activate scrna-integration
  python scripts/make_test_subset.py
"""

import gc
import gzip
import os
import sys
import warnings
from pathlib import Path
from typing import Optional

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

# Source paths (READ-ONLY, must not write)
KIM_DIR = Path("/Users/zhongzishao/Works/GCPL_scRNA/data/Kim_2025_Series GSE295401")
NANCANG_DIR = Path("/Users/zhongzishao/Works/GCPL_scRNA/data/Nancang_2025_GSE249874")
NOWICKI_PATH = "/Users/zhongzishao/Works/GCPL_scRNA/data/Nowicki_Osuch_et_al_2023/Nowicki_Osuch_et_al_2023.h5ad"
TSUBOSAKA_PATH = "/Users/zhongzishao/Works/GCPL_scRNA/data/Tsubosaka_2023_HCA/data_9_9_annotated_seurat_all_ut.rds"
YUE_DIR = Path("/Users/zhongzishao/Works/GCPL_scRNA/data/Yue_SSK_2025_GSE210991")
FIXTURE_B_SOURCE = "/Users/zhongzishao/Works/GCPL_scRNA/results/data_objects/02_qc_filtered_data.h5ad"

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
                        print(f"    Raw: no matching barcodes found (skipping raw)")
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
    all_adatas = []
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
        print("  Fixture A (4/5 sources) + Fixture B ready for stage1-3+ testing.")
    else:
        print("  Some sources failed unexpectedly — check output above.")


if __name__ == "__main__":
    main()
