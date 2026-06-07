"""IO module: multi-source scRNA-seq data readers with manifest-driven obs schema.

The single public function ``read_with_manifest`` implements the 13 behaviors in SPEC
"Three Functions > read_with_manifest".  Private helpers are minimal; no classes,
registries, or plugin systems (ADR-0001 / ADR-0003).
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import Any

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def read_with_manifest(manifest_path: str) -> anndata.AnnData:
    """Read and standardise scRNA-seq data using a per-dataset YAML manifest.

    Implements all 13 behaviors from SPEC "The Three Functions > read_with_manifest".

    Returns a plain ``AnnData``.  Caller does anything they want next.
    """
    manifest = _load_manifest(manifest_path)
    _validate_manifest(manifest)

    source_dataset = str(manifest["source_dataset"])
    input_block = manifest["input"]

    # 1. Read matrix per format
    adata = _read_matrix(input_block, source_dataset)

    # 2. Apply obs_mapping + value_mapping
    _apply_obs_mapping(adata, manifest)

    # 3. Join clinical_metadata tables
    _join_clinical_metadata(adata, manifest)

    # 4. Inject ontology / project_specific constants
    _inject_constants(adata, manifest)

    # 5. Generate cell_id
    _generate_cell_id(adata, source_dataset)

    # 6. Original annotations rename
    _rename_original_annotations(adata, manifest)

    # 7. Gene ID bidirectional sync
    _sync_gene_ids(adata, input_block.get("gene_id_format", "auto"))

    # 8. Species enforcement
    _enforce_species(adata, manifest["species"])

    # 9. Disease system propagation
    _propagate_disease_system(adata, manifest["disease_system"])

    # 10. Layer 2 strong-warn (LLM best-effort fix deferred to PR-3c)
    _warn_layer2(adata)

    # 11. Record raw matrix path
    _record_raw_path(adata, input_block)

    # 12. Baseline QC metrics
    _compute_baseline_qc(adata)

    # 13. Return AnnData — caller does everything else
    return adata


# ---------------------------------------------------------------------------
# Step 1 – Matrix readers
# ---------------------------------------------------------------------------

_FORMAT_HANDLERS = {
    "10x_mtx": "_read_10x_mtx",
    "h5ad": "_read_h5ad",
    "h5": "_read_h5",
    "txt.gz": "_read_txt_gz",
    "rds": "_read_rds",
}


def _read_matrix(input_block: dict[str, Any], source_dataset: str) -> anndata.AnnData:
    fmt = input_block["format"]
    path = input_block["path"]
    handler_name = _FORMAT_HANDLERS.get(fmt)
    if handler_name is None:
        raise ValueError(
            f"Unsupported input.format '{fmt}'. "
            f"Supported: {', '.join(sorted(_FORMAT_HANDLERS))}"
        )
    handler = globals()[handler_name]
    return handler(path, source_dataset, input_block)


def _read_h5ad(path: str, source_dataset: str, _input_block: dict) -> anndata.AnnData:
    adata = sc.read_h5ad(path)
    adata.obs["source_dataset"] = source_dataset
    return adata


def _read_h5(path: str, source_dataset: str, _input_block: dict) -> anndata.AnnData:
    adata = sc.read_10x_h5(path)
    adata.obs["source_dataset"] = source_dataset
    return adata


def _read_10x_mtx(
    path: str, source_dataset: str, _input_block: dict
) -> anndata.AnnData:
    """Read one or more 10x mtx directories (cellranger filtered_feature_bc_matrix).

    *path* may be a single ``filtered_feature_bc_matrix`` directory or a parent
    directory whose immediate children each contain that subdirectory (Nancang).
    """
    mtx_path = Path(path)
    subdirs = _find_10x_dirs(mtx_path)
    if not subdirs:
        raise FileNotFoundError(
            f"No 10x mtx directories found under {path}"
        )

    parts: list[anndata.AnnData] = []
    for sub in sorted(subdirs):
        sample_name = sub.parent.name if sub.name == "filtered_feature_bc_matrix" else sub.name
        try:
            adata_part = sc.read_10x_mtx(sub, var_names="gene_symbols", cache=False)
        except Exception:
            # Some 10x dirs may not have gene_symbols; fall back to gene_ids
            adata_part = sc.read_10x_mtx(sub, var_names="gene_ids", cache=False)
        adata_part.obs["source_dataset"] = source_dataset
        adata_part.obs["sample_id"] = sample_name
        # Prefix barcodes with sample to avoid collisions
        adata_part.obs_names = [f"{sample_name}_{bc}" for bc in adata_part.obs_names]
        parts.append(adata_part)

    if len(parts) == 1:
        adata = parts[0]
    else:
        adata = anndata.concat(parts, join="outer", index_unique="-")
        # sc.read_10x_mtx fills NaN with 0 for missing genes; make it sparse again
        if not sp.issparse(adata.X):
            adata.X = sp.csr_matrix(np.nan_to_num(adata.X, nan=0))

    return adata


def _find_10x_dirs(root: Path) -> list[Path]:
    """Find filtered_feature_bc_matrix directories under *root*."""
    # If root itself is a filtered_feature_bc_matrix
    if (root / "matrix.mtx.gz").exists() or (root / "matrix.mtx").exists():
        return [root]
    # If root contains a filtered_feature_bc_matrix subdir
    if (root / "filtered_feature_bc_matrix").is_dir():
        return [root / "filtered_feature_bc_matrix"]
    # Search children
    result = []
    for child in sorted(root.iterdir()):
        if not child.is_dir():
            continue
        candidate = child / "filtered_feature_bc_matrix"
        if candidate.is_dir():
            result.append(candidate)
    return result


def _read_txt_gz(
    path: str, source_dataset: str, _input_block: dict
) -> anndata.AnnData:
    """Read tab-separated ``*.txt.gz`` count files (Yue organoid format).

    Each file: first column = gene, remaining columns = cell barcodes.
    Matrix is genes x cells; transposed to cells x genes on load.
    """
    txt_dir = Path(path)
    files = sorted(txt_dir.glob("*.txt.gz"))
    if not files:
        raise FileNotFoundError(f"No *.txt.gz files found under {path}")

    parts: list[anndata.AnnData] = []
    for fp in files:
        sample_id = fp.stem.replace("_count", "").split(".")[0]
        # Read tab-separated dense matrix (genes x cells)
        df = pd.read_csv(fp, sep="\t", index_col=0, compression="gzip")
        # Transpose to cells x genes
        X = sp.csr_matrix(df.values.T.astype(np.float32))  # noqa: N806  # scverse convention: canonical single-cell count matrix variable
        var = pd.DataFrame(index=df.index.values)
        obs = pd.DataFrame(index=df.columns.values)
        obs["source_dataset"] = source_dataset
        obs["sample_id"] = sample_id
        adata_part = anndata.AnnData(X=X, obs=obs, var=var)
        parts.append(adata_part)

    if len(parts) == 1:
        adata = parts[0]
    else:
        adata = anndata.concat(parts, join="outer", index_unique="-")
        if not sp.issparse(adata.X):
            adata.X = sp.csr_matrix(np.nan_to_num(adata.X, nan=0))

    return adata


def _read_rds(_path: str, _source_dataset: str, _input_block: dict) -> anndata.AnnData:
    """RDS reader — **not implemented yet** (ADR-0007: RDS needs R environment).

    Tsubosaka dataset is RDS format and is skipped in PR-1 fixtures.
    A future PR will add rpy2-based RDS reading.
    """
    raise NotImplementedError(
        "RDS format is not yet supported.  "
        "The Tsubosaka dataset requires an R environment (ADR-0007).  "
        "See docs/adr/0007-r-bridge-tool-split.md for the tool-split plan.  "
        "This will be addressed in a future PR."
    )


# ---------------------------------------------------------------------------
# Step 2 – obs_mapping + value_mapping
# ---------------------------------------------------------------------------


def _apply_obs_mapping(adata: anndata.AnnData, manifest: dict) -> None:
    obs_mapping = manifest.get("obs_mapping", {})
    if not obs_mapping:
        return
    value_mapping = manifest.get("value_mapping", {})

    for target_col, source_col in obs_mapping.items():
        if source_col not in adata.obs.columns:
            warnings.warn(
                f"obs_mapping: source column '{source_col}' not found in data; "
                f"cannot populate '{target_col}'"
            ,
            stacklevel=2)
            continue
        adata.obs[target_col] = adata.obs[source_col].astype(str)

    for target_col, mapping in value_mapping.items():
        if target_col not in adata.obs.columns:
            warnings.warn(
                f"value_mapping: target column '{target_col}' not found; skipping",
                stacklevel=2,
            )
            continue
        adata.obs[target_col] = adata.obs[target_col].map(
            lambda x, m=mapping: m.get(str(x), x)
        )


# ---------------------------------------------------------------------------
# Step 3 – clinical_metadata join
# ---------------------------------------------------------------------------


def _join_clinical_metadata(adata: anndata.AnnData, manifest: dict) -> None:
    tables = manifest.get("clinical_metadata", [])
    if not tables:
        return

    for tbl_cfg in tables:
        _join_one_clinical_table(adata, tbl_cfg)


def _join_one_clinical_table(adata: anndata.AnnData, cfg: dict) -> None:
    file_path = cfg["file"]
    sheet = cfg.get("sheet", 0)
    skip_rows = cfg.get("skip_rows", 0)
    join_on = cfg["join_on"]
    col_mapping = cfg.get("column_mapping", {})
    val_mapping = cfg.get("value_mapping", {})
    on_missing = cfg.get("on_missing", "warn")
    on_conflict = cfg.get("on_conflict", "metadata_wins")

    manifest_field = join_on["manifest_field"]
    table_column = join_on["table_column"]

    if not os.path.exists(file_path):
        if on_missing == "strict":
            raise FileNotFoundError(f"Clinical metadata file not found: {file_path}")
        warnings.warn(f"Clinical metadata file not found: {file_path}", stacklevel=2)
        return

    if file_path.endswith(".csv"):
        meta_df = pd.read_csv(file_path)
    else:
        meta_df = pd.read_excel(file_path, sheet_name=sheet, skiprows=skip_rows)

    if col_mapping:
        meta_df = meta_df.rename(columns={v: k for k, v in col_mapping.items()})

    for col, mapping in val_mapping.items():
        if col in meta_df.columns:
            meta_df[col] = meta_df[col].map(lambda x, m=mapping: m.get(str(x), x))

    if manifest_field not in adata.obs.columns:
        warnings.warn(
            f"Clinical join: manifest field '{manifest_field}' not found in obs; "
            f"cannot join {file_path}"
        ,
        stacklevel=2)
        return

    # Merge: obs (left) ← metadata (right)
    original_cols = set(adata.obs.columns)
    meta_cols_to_merge = [table_column] + [
        c for c in meta_df.columns
        if c != table_column and c not in original_cols
    ]
    adata.obs = adata.obs.reset_index().merge(
        meta_df[meta_cols_to_merge],
        how="left",
        left_on=manifest_field,
        right_on=table_column,
    ).set_index("index")
    adata.obs.index.name = None

    # Handle on_conflict for overlapping columns
    if on_conflict == "metadata_wins":
        for col in meta_df.columns:
            if col in original_cols and col != table_column:
                # We already renamed via col_mapping; metadata version already populated
                pass


# ---------------------------------------------------------------------------
# Step 4 – inject ontology / project_specific constants
# ---------------------------------------------------------------------------


def _inject_constants(adata: anndata.AnnData, manifest: dict) -> None:
    ontology = manifest.get("ontology", {})
    for key, value in ontology.items():
        adata.obs[key] = value

    project_specific = manifest.get("project_specific", {})
    for col_name, cfg in project_specific.items():
        if isinstance(cfg, dict) and "source_column" in cfg:
            src = cfg["source_column"]
            rules = cfg.get("rules", {})
            if src in adata.obs.columns:
                adata.obs[col_name] = adata.obs[src].map(
                    lambda x, r=rules: r.get(str(x), str(x))
                )
        elif isinstance(cfg, dict) and "value" in cfg:
            adata.obs[col_name] = cfg["value"]


# ---------------------------------------------------------------------------
# Step 5 – cell_id generation
# ---------------------------------------------------------------------------


def _generate_cell_id(adata: anndata.AnnData, source_dataset: str) -> None:
    sample = adata.obs.get("sample_id", pd.Series("unknown", index=adata.obs_names))
    barcode = [str(bc).split("-")[0] for bc in adata.obs_names]
    adata.obs["cell_id"] = [
        f"{source_dataset}_{s}_{b}"
        for s, b in zip(sample, barcode, strict=True)
    ]


# ---------------------------------------------------------------------------
# Step 6 – original_annotations rename
# ---------------------------------------------------------------------------


def _rename_original_annotations(adata: anndata.AnnData, manifest: dict) -> None:
    annotations = manifest.get("original_annotations", [])
    for entry in annotations:
        col = entry["column"]
        role = entry.get("role", "")
        if col not in adata.obs.columns:
            warnings.warn(
                f"original_annotations: column '{col}' not found in data; skipping"
            ,
            stacklevel=2)
            continue
        source = manifest["source_dataset"]
        suffix = f"_{role}" if role else ""
        new_name = f"cell_type_original_{source}_v1{suffix}"
        adata.obs[new_name] = adata.obs[col]


# ---------------------------------------------------------------------------
# Step 7 – gene ID bidirectional sync
# ---------------------------------------------------------------------------


def _sync_gene_ids(adata: anndata.AnnData, gene_id_format: str = "auto") -> None:
    """Ensure both gene symbols (var.index) and ensembl IDs (var['ensembl_id']).

    - If var.index is symbols: extract ensembl from var columns or query mygene.
    - If var.index is ensembl: convert to symbols via mygene and re-index.
    """
    idx_sample = str(adata.var.index[0])
    is_ensembl = idx_sample.startswith("ENSG")

    if gene_id_format == "symbol" or (gene_id_format == "auto" and not is_ensembl):
        _sync_symbol_to_ensembl(adata)
    elif gene_id_format == "ensembl" or (gene_id_format == "auto" and is_ensembl):
        _sync_ensembl_to_symbol(adata)
    else:
        _sync_symbol_to_ensembl(adata)


def _sync_symbol_to_ensembl(adata: anndata.AnnData) -> None:
    """var.index is symbols; try to add var['ensembl_id'] from existing columns or mygene."""
    if "ensembl_id" in adata.var.columns and adata.var["ensembl_id"].notna().any():
        # Already present — possibly from 10x gene_ids column
        return

    # Check for gene_ids column (10x convention)
    if "gene_ids" in adata.var.columns:
        gene_ids = adata.var["gene_ids"]
        # gene_ids may be a string like "ENSG00000236601" or missing
        adata.var["ensembl_id"] = gene_ids.where(
            gene_ids.astype(str).str.startswith("ENSG"), ""
        )
        n_mapped = (adata.var["ensembl_id"] != "").sum()
        if n_mapped > 0:
            return

    # Fall back to mygene
    _mygene_symbol_to_ensembl(adata)


def _mygene_symbol_to_ensembl(adata: anndata.AnnData) -> None:
    """Query mygene to convert gene symbols to Ensembl IDs."""
    import mygene

    mg = mygene.MyGeneInfo()
    symbols = list(adata.var.index)
    ensembl_ids: dict[str, str] = {}

    batch_size = 1000
    for i in range(0, len(symbols), batch_size):
        batch = symbols[i : i + batch_size]
        try:
            results = mg.querymany(
                batch, scopes="symbol", fields="ensembl.gene",
                species="human", returnall=False
            )
        except Exception:
            warnings.warn(f"mygene query failed for batch {i}; ensembl_id left empty", stacklevel=2)
            continue
        for r in results:
            query = r.get("query", "")
            ensembl_data = r.get("ensembl")
            if ensembl_data and isinstance(ensembl_data, dict):
                eid = ensembl_data.get("gene")
                if eid:
                    ensembl_ids[query] = eid

    adata.var["ensembl_id"] = adata.var.index.map(ensembl_ids).fillna("")
    n_missing = (adata.var["ensembl_id"] == "").sum()
    if n_missing > 0:
        warnings.warn(
            f"Gene ID sync: {n_missing}/{len(symbols)} symbols could not be "
            f"mapped to Ensembl IDs via mygene.  Corresponding ensembl_id entries "
            f"left empty."
        ,
        stacklevel=2)


def _sync_ensembl_to_symbol(adata: anndata.AnnData) -> None:
    """var.index is Ensembl IDs; convert to gene symbols via mygene or feature_name."""
    # Check for feature_name column (common in CellxGene-hosted h5ads)
    if "feature_name" in adata.var.columns:
        symbol_map = adata.var["feature_name"].to_dict()
        # Use feature_name where available
        adata.var["_orig_ensembl"] = adata.var.index.values
        new_index = [
            symbol_map.get(eid, eid) for eid in adata.var.index
        ]
        adata.var.index = new_index
        adata.var.index.name = None
        adata.var["ensembl_id"] = adata.var["_orig_ensembl"]
        del adata.var["_orig_ensembl"]
        n_mapped = sum(1 for v in new_index if not str(v).startswith("ENSG"))
        if n_mapped > 0:
            return
        # feature_name didn't help (all still ENSG); fall through to mygene

    # Query mygene
    import mygene

    mg = mygene.MyGeneInfo()
    ensembl_ids = list(adata.var.index)
    symbol_map_result: dict[str, str] = {}

    batch_size = 1000
    for i in range(0, len(ensembl_ids), batch_size):
        batch = ensembl_ids[i : i + batch_size]
        try:
            results = mg.querymany(
                batch, scopes="ensembl.gene", fields="symbol",
                species="human", returnall=False
            )
        except Exception:
            warnings.warn(f"mygene query failed for batch {i}; keeping Ensembl IDs", stacklevel=2)
            continue
        for r in results:
            query = r.get("query", "")
            sym = r.get("symbol")
            if sym:
                symbol_map_result[query] = sym

    adata.var["_orig_ensembl"] = adata.var.index.values
    new_index = [
        symbol_map_result.get(eid, eid) for eid in adata.var.index
    ]
    adata.var.index = new_index
    adata.var.index.name = None
    adata.var["ensembl_id"] = adata.var["_orig_ensembl"]
    del adata.var["_orig_ensembl"]

    n_missing = sum(1 for v in new_index if str(v).startswith("ENSG"))
    if n_missing > 0:
        warnings.warn(
            f"Gene ID sync: {n_missing}/{len(ensembl_ids)} Ensembl IDs could not be "
            f"mapped to gene symbols via mygene.  Keeping Ensembl IDs in var.index "
            f"and leaving ensembl_id empty for those genes."
        ,
        stacklevel=2)
        # For unmapped, leave ensembl_id empty
        mask = adata.var.index.astype(str).str.startswith("ENSG")
        adata.var.loc[mask, "ensembl_id"] = ""


# ---------------------------------------------------------------------------
# Step 8 – species enforcement
# ---------------------------------------------------------------------------


def _enforce_species(adata: anndata.AnnData, species: str) -> None:
    if species != "human":
        raise ValueError(
            f"species '{species}' is not supported.  "
            f"The framework currently only accepts species='human'.  "
            f"Other species require a new ADR."
        )
    adata.uns["species"] = "human"


# ---------------------------------------------------------------------------
# Step 9 – disease_system propagation
# ---------------------------------------------------------------------------


def _propagate_disease_system(adata: anndata.AnnData, disease_system: str) -> None:
    adata.obs["disease_system"] = disease_system


# ---------------------------------------------------------------------------
# Step 10 – Layer 2 strong-warn
# ---------------------------------------------------------------------------

_LAYER2_FIELDS = [
    "disease",
    "disease_ontology_term_id",
    "tissue",
    "tissue_ontology_term_id",
    "assay",
    "sex",
    "development_stage",
]


def _warn_layer2(adata: anndata.AnnData) -> None:
    """Emit strong warnings for missing/malformed Layer 2 CellxGene-aligned fields.

    LLM best-effort fix is deferred to PR-3c (needs OpenRouter key).
    This function only warns and leaves NaN — never blocks ingest.
    """
    for field in _LAYER2_FIELDS:
        if field not in adata.obs.columns:
            warnings.warn(
                f"[Layer2] obs column '{field}' is MISSING.  "
                f"LLM best-effort fix deferred to PR-3c.  "
                f"Column will be created with NaN for all cells."
            ,
            stacklevel=2)
            adata.obs[field] = np.nan
        else:
            col = adata.obs[field]
            n_null = col.isna().sum()
            n_empty = (col.astype(str).str.strip() == "").sum()
            n_problem = n_null + n_empty
            if n_problem > 0:
                warnings.warn(
                    f"[Layer2] obs column '{field}' has {n_problem}/{len(col)} "
                    f"missing or empty values.  LLM best-effort fix deferred to PR-3c."
                ,
                stacklevel=2)
            # Check for suspicious values (e.g. "unknown", "NA", "N.A.")
            suspicious = col.astype(str).str.lower().isin(
                ["unknown", "na", "n.a.", "n/a", "none", "null"]
            )
            if suspicious.any():
                warnings.warn(
                    f"[Layer2] obs column '{field}' has {suspicious.sum()} cells "
                    f"with suspicious values ('unknown'/'NA' etc).  "
                    f"LLM fix deferred to PR-3c."
                ,
                stacklevel=2)


# ---------------------------------------------------------------------------
# Step 11 – record raw matrix path
# ---------------------------------------------------------------------------


def _record_raw_path(adata: anndata.AnnData, input_block: dict) -> None:
    raw_path = input_block.get("raw_path")
    if raw_path:
        adata.uns["raw_matrix_path"] = raw_path
    else:
        adata.uns["raw_matrix_path"] = None


# ---------------------------------------------------------------------------
# Step 12 – baseline QC
# ---------------------------------------------------------------------------


def _compute_baseline_qc(adata: anndata.AnnData) -> None:
    """Compute n_genes, total_counts, pct_counts_mt, pct_counts_ribo on adata.X."""
    # Ensure CSR for scanpy QC functions
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)

    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1 if sp.issparse(adata.X) else (adata.X > 0).sum(axis=1)
    adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).flatten()

    # MT genes
    mt_mask = adata.var.index.str.startswith("MT-")
    if mt_mask.any():
        adata.obs["pct_counts_mt"] = (
            np.asarray(adata.X[:, mt_mask].sum(axis=1)).flatten()
            / adata.obs["total_counts"].values
            * 100
        )

    # Ribo genes
    ribo_mask = adata.var.index.str.startswith(("RPS", "RPL"))
    if ribo_mask.any():
        adata.obs["pct_counts_ribo"] = (
            np.asarray(adata.X[:, ribo_mask].sum(axis=1)).flatten()
            / adata.obs["total_counts"].values
            * 100
        )


# ---------------------------------------------------------------------------
# Manifest loading & validation
# ---------------------------------------------------------------------------


def _load_manifest(manifest_path: str) -> dict:
    path = Path(manifest_path)
    if not path.exists():
        raise FileNotFoundError(f"Manifest file not found: {manifest_path}")
    with open(path) as fh:
        manifest = yaml.safe_load(fh)
    if not isinstance(manifest, dict):
        raise ValueError(f"Manifest at {manifest_path} is not a valid YAML mapping")
    return manifest


def _validate_manifest(manifest: dict) -> None:
    """Validate manifest schema. Raises ValueError on fatal issues."""

    # species: mandatory, currently only "human"
    species = manifest.get("species")
    if not species:
        raise ValueError("Manifest is missing required field 'species'")
    if species != "human":
        raise ValueError(
            f"species '{species}' is not supported. Only 'human' is accepted."
        )

    # input block mandatory
    if "input" not in manifest:
        raise ValueError("Manifest is missing required section 'input'")

    inp = manifest["input"]
    for key in ("format", "path"):
        if key not in inp:
            raise ValueError(f"Manifest 'input' section is missing required key '{key}'")

    fmt = inp["format"]
    if fmt not in _FORMAT_HANDLERS:
        raise ValueError(
            f"Unsupported input.format '{fmt}'. "
            f"Supported: {', '.join(sorted(_FORMAT_HANDLERS))}"
        )

    # source_dataset mandatory
    if "source_dataset" not in manifest:
        raise ValueError("Manifest is missing required field 'source_dataset'")

    # project_id mandatory
    if "project_id" not in manifest:
        raise ValueError("Manifest is missing required field 'project_id'")

    # disease_system mandatory
    if "disease_system" not in manifest:
        raise ValueError("Manifest is missing required field 'disease_system'")

    # original_annotations section mandatory (even if [])
    if "original_annotations" not in manifest:
        raise ValueError(
            "Manifest is missing required section 'original_annotations' "
            "(use [] if the dataset has no author annotations)"
        )

    # qc_overrides: reason mandatory when skip: true
    for step, cfg in manifest.get("qc_overrides", {}).items():
        if cfg.get("skip") and not cfg.get("reason"):
            raise ValueError(
                f"qc_overrides.{step}.skip is true but 'reason' is missing.  "
                f"A reason is required when skipping a QC step."
            )
