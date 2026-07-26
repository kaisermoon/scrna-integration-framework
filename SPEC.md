# Implementation Specification

Implementation details for the scRNA-seq Integration Framework. Companion to `CONTEXT.md` (project glossary) and `docs/adr/` (architectural decisions). This document specifies *what is built and how it behaves*; CONTEXT.md defines *what terms mean*.

---

## Positioning

**The framework is a tool for biological discovery, not a methods contribution.** Its purpose is to enable PI and students to produce publishable biological findings across disease projects (gastric precancerous lesions, RA, PCOS, VVC, ...) — not to publish a software paper or release a community library. This framing rules out the design pressures that come with public-tool ambitions (broad API stability, exhaustive docs, generic extensibility).

**01 (current)**: Personal research infrastructure. Primary users are PI and supervised students. API may iterate aggressively; no public stability guarantees; lean docs. Function signatures may change between PRs during this stage — each signature change is logged in `_memory.md` so future readers can trace when an API shifted.

**02 (later)**: Re-evaluated when the framework has produced **at least one publishable biological finding from the GCPL pilot** (PR-4 complete + a 07 downstream module completes a real analysis). What "02" means at that point is decided then — it is not committed in advance to be a software paper, a public release, or any specific external deliverable.

### Engineering discipline from day one

These apply from the first PR, regardless of audience size:

- **Modularity** — the two-function framework boundary; ADR-0001 / ADR-0003 / ADR-0004 enforced by code review.
- **Tests** — unit tests on framework code; smoke tests on stage notebooks.
- **Pre-commit** — hooks installed; gitleaks; basic checks; ADR-mandated patterns flagged at commit.
- **Branch protection** — main protected, PR-required, CI-required (already configured).
- **ruff** — Python linter for unused imports / dead code / outdated syntax. Already configured in `pyproject.toml`; runs in CI and pre-commit. Zero learning cost.
- **mypy is NOT used** — scientific Python ecosystem (scanpy / anndata / scvi-tools) has too many polymorphic return types; mypy would generate false positives without catching real bugs. Type annotations are optional, added only when they help readability.

External commitments (full docs site, semver, deprecation policy, PyPI release) are explicitly out of scope. They are not a "later" promise — they are simply not goals.

---

## Architectural Stance: Two Functions + Scorers, Conventions

**The framework code surface is intentionally tiny — designed for non-CS PI/students to read line-by-line (ADR-0009).** Two Python functions cover the genuine gaps in scanpy/anndata/scverse; the scorers module provides directly-callable metric functions (no callbacks, no `sweep()` abstraction); everything else is convention plus reference data, used directly with `pandas` / `pyyaml` / `scanpy` / `anndata` native APIs.

```python
from scrna_integration import read_with_manifest, load_markers
# scorers are imported directly from the module — call in plain for loops:
from scrna_integration.scorers import integration_metrics, clustering_metrics
```

Notebooks import only what they need. Anything else (lineage, disease ontology lookup, QC skip records, run metadata, stage reports) is plain `adata.uns[...]` writes, plain `yaml.safe_load`, plain `sc.pl.*`, or notebook cells PI / students edit directly.

This stance follows from:

- **ADR-0001** — thin framework over scanpy. Don't wrap functions PI already knows.
- **ADR-0003** — plain code over plugin systems. Pass plain functions / dicts, don't build registries.
- **ADR-0004** — framework deletion log. Records the second-round PI review of `references/legacy-GCPL/` showing every notebook is scanpy-native end-to-end.
- **ADR-0005** — `load_markers` as a legitimate third function. Records the ADR-0004 escape-hatch criteria being met for marker-library loading.

### What the framework does NOT provide

The following were considered and rejected as wrapper noise:

| Rejected wrapper | Use this instead |
|---|---|
| `si.tracking.register_run(...)` | `adata.uns["my_run_2026"] = {"params": ..., "summary": ...}` |
| `si.lineage.promote(...)` / `show_lineage` / `show_dependents` | File-naming convention + `adata.uns["status"] = "promoted"` + `glob` over `results/` when needed |
| `si.disease.ancestors(...)` | `yaml.safe_load(open("references/disease_ontology/gastric.yaml"))` then PI walks the dict |
| `si.report.stage*(adata, ...)` | Stage notebooks with `sc.pl.*` + `plt.savefig` cells PI edits directly |
| `si.qc.show_skipped(adata)` | `adata.uns["qc_skipped"]` direct read |
| `si.memory.check(adata)` | code reviewer agent enforces conventions; not a runtime function |
| `si.scorers.qc_balance` | Plain functions in `scrna_integration/scorers.py`, importable when desired but not framework-namespaced |

The framework's `__init__.py` re-exports `read_with_manifest` and `load_markers`. Scorers are imported directly from `scrna_integration.scorers` when needed in notebooks. There is no `si.*` sub-namespace tree.

### "Two functions + scorers" is the current reality, not an eternal promise

The current framework surface reflects today's understanding of where scanpy genuinely leaves a gap. A third (or further) framework function may be added when **all three** of the following are true:

1. **The same boilerplate appears in ≥3 stage notebooks** with non-trivial duplication (not just a 1-liner).
2. **The naive approach has been tried in real PR work and demonstrably caused maintenance burden** — copy-paste drift across notebooks, missed updates when the pattern needs to evolve, etc. "I anticipate this will become annoying" is not enough; the maintenance cost must already be observed (or — as with `load_markers` per ADR-0005 — PI confirms from prior real-world experience that the duplication is observed in lived practice, not anticipated).
3. **A new ADR is written** explaining why this case is an exception to ADR-0001 / 0003 / 0004 — what the naive approach was, why it failed, and what the new function does.

This bar is intentionally high. Most candidate "abstractions" that emerge during PR work fail criterion 2 — perceived future pain is not a license to add framework code. When in doubt, accept the duplication and revisit later.

If the framework grows from 2 functions to 5+, that is a signal the bar is being abused, not a signal the bar is wrong — tighten the bar (e.g. require >=5 notebooks or PI sign-off in addition to the ADR).

---

## Two Functions + Scorers

### `read_with_manifest(manifest_path: str) -> AnnData`

The genuine IO gap: scanpy's readers don't unify cross-source obs schemas, don't join clinical metadata tables, don't apply value mappings. This function:

1. Reads the raw matrix per `input.format` (10x_mtx / h5ad / h5 / rds via rpy2).
2. Applies `obs_mapping` + `value_mapping` from the manifest.
3. Joins each `clinical_metadata` table with the declared `join_on` / `on_missing` / `on_conflict`.
4. Injects ontology and project-specific constants.
5. Auto-generates `cell_id = {source_dataset}_{sample_id}_{barcode}` to prevent barcode collisions.
6. Renames declared `original_annotations` columns to `cell_type_original_{source_dataset}_v1[_{role}]`.
7. **Bidirectional gene ID sync**: detects whether the input `var.index` is gene symbols or Ensembl IDs and ensures both are present after read. `var.index` is canonicalised to gene symbols (matching scanpy convention and `legacy-GCPL`); the corresponding Ensembl IDs are kept as `var["ensembl_id"]`. If the input only has Ensembl IDs, framework converts to symbols (via `mygene` or equivalent) and uses symbols as `var.index`. Untranslatable rows keep an empty `ensembl_id` and emit a warning.
8. **Species enforcement**: requires `species` field in manifest (currently only `"human"` accepted; other species fail loudly). Writes `adata.uns["species"] = "human"` so downstream stages can verify before loading species-specific markers / references.
9. **Disease system propagation**: copies the manifest's top-level `disease_system` field to `obs.disease_system` (Layer 1 required field) for every cell. Cross-project integration analyses across PI's portfolio (gastric / RA / PCOS / VVC) rely on this column.
10. **Layer 2 strong-warn + LLM best-effort fix**: for each of the seven CellxGene-aligned fields (`disease`, `disease_ontology_term_id`, `tissue`, `tissue_ontology_term_id`, `assay`, `sex`, `development_stage`) that is missing or malformed, emits a strong warning and invokes an LLM-assisted disambiguator (looking at obs head + manifest + clinical-metadata head) to propose a value. PI confirms in-notebook; confirmed fixes are written back into `manifest.yaml` for reproducibility.
11. **Records (does not load) raw matrix path** for SoupX: if manifest provides `input.raw_path`, writes it to `adata.uns["raw_matrix_path"]`. 02 SoupX reads it on demand to avoid doubling 01 memory.
12. Computes baseline QC metrics (`n_genes`, `total_counts`, `pct_counts_mt`, `pct_counts_ribo`) on `adata.X` so they are present and column-aligned across all source datasets regardless of preprocessing state. These are 02 prerequisites; computing them at ingest avoids a missing-column branch at 02.
13. Returns plain `AnnData`. **Caller does anything they want next** — `adata.write_h5ad(...)` to checkpoint, scanpy ops to continue, etc.

The function is one file (`src/scrna_integration/io.py`) + minimal helpers. It does **not** call `validate_obs`, does **not** push metadata into a hidden namespace beyond `species` / `raw_matrix_path`, does **not** assign a "stage" tag. PI inspects obs after reading and decides whether the schema is OK.

### Scorers: directly-callable metric functions (no `sweep()` abstraction)

The `scrna_integration.scorers` module provides metric functions with simple signatures
suitable for calling directly in notebook for loops — no callbacks, no `fn=`/`scorer=` parameters.
Designed to be transparently readable by non-CS students (ADR-0009).

```python
from scrna_integration.scorers import integration_metrics, clustering_metrics
```

**Explicit for-loop pattern** (replaces the removed `sweep()` function — ADR-0009):
```python
# 04: iterate embeddings, compute metrics directly
results = []
for rep in use_reps:
    adata_copy = adata.copy()
    sc.pp.neighbors(adata_copy, use_rep=rep)
    sc.tl.umap(adata_copy)
    m = integration_metrics(adata_copy)   # direct call, no callback
    results.append({"use_rep": rep, **m})
import pandas as pd; sweep_df = pd.DataFrame(results)
```

05 uses the same pattern with `clustering_metrics` across resolutions.

Available scorer functions:
- `integration_metrics(adata, batch_key="batch", label_key=None)` — silhouette by batch + celltype
- `clustering_metrics(adata, cluster_key=None, label_key=None)` — silhouette + ARI
- `qc_balance(adata_after, adata_before)` — cell/gene retention after QC
- `annotation_concordance(adata, label_a=None, label_b=None)` — Cohen's kappa

All gracefully handle missing columns/embeddings by returning `_note": NaN`.
scib-metrics is optional — `scib_available` = 0.0 when not installed.

### `load_markers(csv_path, roles=("canonical", "optional")) -> dict`

A marker-library loader for the `references/markers/*.csv` corpus. Justified per ADR-0005: PI's accumulated experience confirms the load + filter + groupby boilerplate appears in 06 annotation, per-cluster profiling, and downstream gene-set scoring notebooks; centralising the `role` semantics (`canonical` / `optional` / `negative`) prevents misuse.

```python
def load_markers(
    csv_path: str,                                      # any path; relative or absolute
    roles: tuple[str, ...] | None = ("canonical", "optional"),
) -> dict[str, list[str]] | dict[str, dict[str, list[str]]]:
    """
    Load a marker CSV and return cell_type → markers, filtered by role.

    Default (canonical + optional) is the gene-set scoring / dotplot common case.
    roles=("negative",)  → negative markers for reverse validation
    roles=None           → full 3-layer dict {cell_type: {canonical: ..., optional: ..., negative: ...}}
    """
```

Typical usage:

```python
markers = load_markers("references/markers/gastric_epithelial.csv")
# → {"SPEM": ["TFF2", "MUC6"], "pit_cell": ["MUC5AC", "TFF1"], ...}
sc.pl.dotplot(adata, var_names=markers, groupby="leiden")
```

The function is one file (`src/scrna_integration/markers.py`). It accepts any path so PI can keep ad-hoc marker files outside `references/markers/` if needed for one-off analyses.

---

## Pipeline Stages — File-naming Convention Only

The pipeline is a **file-naming convention**, not a framework feature. There is no `Stage` class, no stage entry function, no banner output, no automatic upstream selection.

```
stage 0    raw/                                          external raw data (cellranger / h5ad / RData), read-only
01    results/01_{dataset}_v{N}.h5ad            per-dataset independent QC (MAD adaptive thresholds / scrublet / SoupX / cell cycle / complexity / special gene markers)
02    results/02_merged_v{N}.h5ad              explicit anndata.concat(join="inner") + cross-dataset diagnostics
03    results/03_normalized_v{N}.h5ad          normalize + log1p + batch-aware HVG + HVG exclusion list (+ optional Pearson residuals)
04    results/04_embedded_v{N}.h5ad            PCA / Harmony / scVI / scANVI (elbow plot + N_NEIGHBORS sweep + HARMONY_THETA sweep + integration metrics)
05    results/05_clustered_v{N}.h5ad           multi-resolution Leiden + clustering stability (subsample ARI) + marker gene preview
06    results/06_annotated_v{N}.h5ad           multi-method annotation + cross-method comparison
06c  results/06c_subset_v{N}.h5ad              subset re-cluster (re-run 03-5 on a subset; e.g. all T cells)
07    results/downstream/{module}_v{N}.h5ad    downstream modules (one h5ad per module, see below)
```

### Per-dataset QC (Stage 01)

Each source dataset gets its own notebook under `notebooks/01_per_dataset/`. A dataset notebook independently runs `read_with_manifest` then performs the complete QC pipeline (MAD-based adaptive thresholds, scrublet doublet detection, SoupX ambient correction, cell cycle scoring, gene complexity, special gene marker profiling) before writing a per-dataset checkpoint h5ad. This design (ADR-0011) replaces the former monolithic `01_loaded` + `02_qcd` approach, enabling:

- **Dataset-specific QC thresholds** — tissue biopsies vs organoids have fundamentally different MT% baselines; each gets independent adaptive thresholds
- **Inspect-before-merge** — researchers can examine each dataset's quality distribution independently before deciding to merge
- **Clear QC strategy per dataset** — recorded in `adata.uns["qc_report_v1"]` with `strategy` field (adaptive / fixed / skip), aggregated at merge time

### Merge (Stage 02)

`02_merged.ipynb` loads all per-dataset h5ads and performs explicit `anndata.concat(join="inner", label="source_dataset", keys=...)`. It validates cell_id uniqueness, runs cross-dataset QC diagnostics (violin plots and summary tables per source_dataset), and aggregates per-dataset `uns["qc_report_v1"]` records into a unified `uns["merge_report_v1"]`. The merged AnnData serves as the single upstream for stages 03-15.

### Scalar-or-Sweep: Dual-mode Parameter Design (Stages 03-05)

Stages 03 through 05 support a **dual-mode parameter pattern**: write a single scalar value and the notebook executes that value directly; write a Python list and the notebook automatically sweeps over all values, producing comparison visualisations. There is no `sweep()` function or framework abstraction — each stage notebook uses plain `isinstance(param, list)` branching with explicit for loops (ADR-0009).

**Design goals**:
- **Single value mode** — the happy path. One cell, one parameter, one result. Students reading the notebook see one clean execution path.
- **List/sweep mode** — exploratory. The same notebook cell runs a for loop over multiple values, computes comparison metrics (Jaccard overlap, UMAP grids, silhouette scores, etc.), and produces visualisations to help PI choose the best value.
- **No framework magic** — the branching is visible Python in the notebook: `_values = N_TOP_GENES if isinstance(N_TOP_GENES, list) else [N_TOP_GENES]`. A non-CS student reading the code sees exactly how the sweep works.

**Parameters supporting Scalar-or-Sweep**:

| Stage | Parameter | Single example | Sweep example | Comparison output |
|-------|-----------|---------------|---------------|-------------------|
| 03 | `N_TOP_GENES` | `2000` | `[1500, 2000, 3000, 4000]` | Jaccard overlap heatmap of HVG sets |
| 03 | `HVG_FLAVOR` | `"seurat"` | `["seurat", "seurat_v3"]` | Jaccard overlap heatmap (crossed with `N_TOP_GENES`) |
| 04 | `N_NEIGHBORS` | `15` | `[10, 15, 20, 30]` | UMAP grid coloured by batch_key, comparing local-vs-global structure |
| 04 | `HARMONY_THETA` | `2.0` | `[1, 2, 3]` | UMAP grid comparing batch correction strength (mixing vs over-correction) |
| 05 | `RESOLUTIONS` | `[0.2, 0.4, ..., 2.0]` | (always a list) | Cluster count vs resolution plots, silhouette scores, cluster size distributions |

**Mode switching**: to switch from single-value to sweep mode, replace `N_TOP_GENES = 2000` with `N_TOP_GENES = [1500, 2000, 3000]` in the PARAMS cell. No other code changes needed. The notebook's branching logic handles both paths.

**Sweep comparison outputs**:
- All sweep results are saved under `results/figures/sweep_{stage}/` (e.g. `sweep_04/` for embedding comparisons).
- Each sweep writes a markdown report summarising the comparison table and key findings.
- Sweep outputs are diagnostic aids — they inform PI's parameter choice but are not themselves pipeline checkpoints. The chosen parameter value goes into the next stage's notebook PARAMS cell.

### 04 sweep includes integration QC

No embedding method is best on every dataset, so the workflow is to **run all candidates, then compare** — never to pick a default method up front. Candidates are peers: `X_pca` (baseline), `X_pca_harmony`, `X_scVI`, `X_scANVI` (when a labelled reference exists), cellxgene_census pretrained scVI (optional), and any future method PI adds (e.g. scPoli / scArches-family methods). Each writes its own `obsm` slot; none is privileged.

Two decision inputs, both first-class:

1. **Visual inspection (primary in PI's actual workflow)** — `sc.pl.umap` for every candidate embedding, coloured by `sample_id`, `batch`, and `cell_type` (original-author or marker-based). PI eyeballs how well each integrates batches without over-mixing biologically distinct cells. This visual read is the main basis for the decision, not a formality.
2. **Integration metrics** — batch entropy / iLISI / scIB suite (silhouette by batch / silhouette by celltype / kBET / graph connectivity / etc.), embedded in the 04 sweep report via `sweep(fn=embedding_method, ..., scorer=integration_metrics)`. The metric table corroborates or complicates the visual read.

PI weighs both and marks one embedding `promoted`. cellxgene_census pretrained scVI and scANVI are optional candidate cells (each with its prerequisite noted — census model coverage for the tissue / a labelled reference atlas respectively), not framework-prescribed steps.

**Adding an embedding method (extension pattern).** Same scanpy-native "parallel slots" mechanism as 05 clustering: a new method writes `adata.obsm["X_{method}"]` and is appended to the 04 notebook's `use_reps` list in the explicit for loop — one added cell, no framework change. This is how PI plugs in new integration methods as they are published.

### 06c subset analysis

After 06 PI may want to refine a specific cell-type subset (all T cells, all epithelial cells, all SPEM-spectrum cells, ...). 06c re-runs 03 (HVG re-selection on the subset) → 04 (re-embedding) → 05 (re-clustering) → 06 (re-annotation) on the subset. Output goes to `06c_subset_v{N}.h5ad` with `adata.uns["subset_of"] = "results/gcpl_06_annotated_v1.h5ad"` and `adata.uns["subset_filter"] = "cell_type_final_v1.isin(['CD4 T', 'CD8 T', 'Treg'])"` for traceability.

06c is implemented as a single notebook (`06c_subset.ipynb`) that calls the same scanpy / framework API as stages 3-6 on the subset. Re-clustering parameters are project-specific (the T-cell subset usually needs different HVG / resolution than the global UMAP).

### 07 downstream modules

07 is a fan-out — multiple downstream modules consume 06 (or 06c) output independently and run in parallel. Each module has its own h5ad and its own notebook:

| Module | Notebook | Tool stack | Status | Student-code reference (re-implement per ADR-0008) |
|---|---|---|---|---|
| Differential expression (per-cluster, scanpy native) | `notebooks/07_downstream/deg.ipynb` | `sc.tl.rank_genes_groups` | **PR-3 in scope** | legacy-GCPL `08_differential_expression.ipynb` |
| Pseudobulk DEG (cross-condition: disease vs control) | `notebooks/07_downstream/pseudobulk_deg.ipynb` | DESeq2 (subprocess Rscript) | **PR-3 in scope** | — (decoupler pseudobulk + DESeq2 contrast) |
| CNV inference | `notebooks/07_downstream/cnv.ipynb` | infercnvpy (pure Python) | **PR-3 in scope** | `workflow_for_pseudotime/4.2_*` (gene positions) + `4.3_*` (infercnvpy run) |
| Pathway enrichment | `notebooks/07_downstream/pathway.ipynb` | GSEApy / decoupler / Reactome | PR-5+ | CD4 deep-analysis template GSEA part (reference only) |
| Pseudotime + root identification | `notebooks/07_downstream/pseudotime.ipynb` | CytoTRACE (cellrank) / Monocle3 (Rscript) / transcriptome entropy (numpy) | PR-5+ | `4.3_*` (entropy/CytoTRACE/Monocle3 export) + `4.4_*`/`4.5_*` (multi-metric root id) + `11.2_*` (lineage Monocle3) |
| GRN | `notebooks/07_downstream/grn.ipynb` | pySCENIC | PR-5+ | — |
| Cell communication | `notebooks/07_downstream/cell_communication.ipynb` | CellChat / NicheNet / CellPhoneDB | PR-5+ | CD4 deep-analysis LR part (reference only) |
| Differential abundance | `notebooks/07_downstream/abundance.ipynb` | scCODA / Milo | PR-5+ | `11_all_celltype_proportion_analyse.ipynb` (scCODA + Mann-Whitney + Cliff's delta + effect sizes) |
| Gene co-expression modules | `notebooks/07_downstream/gene_modules.ipynb` | hdWGCNA (subprocess Rscript) | PR-5+ | — |

The 3 modules in PR-3 scope (DEG / pseudobulk DEG / CNV) form the **minimum set needed to reach the GCPL pilot's first biological finding** — distinguishing tumor vs normal cells, finding genes differential along the CAG → IM → dysplasia axis, and identifying disease-vs-control pseudobulk DEGs. The remaining 6 modules are added in subsequent PRs as the GCPL analysis surfaces specific scientific questions.

Multiple method runs at any stage coexist in the same h5ad as additional `obsm` / `obs` / `var` slots, following scanpy's natural anndata usage. Re-runs bump the version suffix; old versions remain on disk for rollback.

### Lineage = file naming + uns annotation

A re-run writes its own provenance directly into `adata.uns`:

```python
adata.uns["status"] = "promoted"          # or "experimental" / "deprecated"
adata.uns["upstream"] = ["results/gcpl_03_normalized_v1.h5ad"]
adata.uns["notes"] = "tried scVI in place of Harmony; not promoted"
```

PI inspects these dict entries directly. There is no API for `promote()` / `deprecate()` / `show_dependents()`. When PI wants to know "what depends on 04_v1", they `glob` over `results/` and read each file's `adata.uns["upstream"]`.

### Default upstream selection

There is no automatic default. Every notebook has a parameters cell where the upstream path is **explicitly named**:

```python
# === PARAMS ===
upstream_path = "results/gcpl_03_normalized_v1.h5ad"
output_version = 2
```

This was originally proposed as a framework banner — replaced by the simpler convention that every stage notebook starts with a `# === PARAMS ===` cell so the upstream is visible at the top of the file. No magic, no resolution rule.

---

## Run Metadata — Plain `adata.uns` Writes

PI / agents writing notebook code record run metadata directly. There is **no framework function for this** — `adata.uns[key] = {...}` is the entire convention.

### Three concrete examples from real workflow

**Example 1 — recording a QC filter run**:

```python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=6000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()

adata.uns["filter_v1"] = {
    "params":  {"min_genes": 200, "max_genes": 6000, "max_pct_mt": 20},
    "cells_in":  50000,
    "cells_out": 42000,
    "timestamp": pd.Timestamp.now().isoformat(),
}
```

**Example 2 — recording a Harmony embedding run** at 04:

```python
sce.pp.harmony_integrate(adata, key="batch")

adata.uns["harmony_v1"] = {
    "method":    "harmony",
    "batch_key": "batch",
    "n_pcs":     30,
    "obsm_key":  "X_pca_harmony",   # cross-reference where harmony's output landed
    "timestamp": pd.Timestamp.now().isoformat(),
}
```

**Example 3 — recording PI's `cell_type_final` decision** at 06:

```python
adata.obs["cell_type_final_v1"] = adata.obs["leiden_res_1.0"].map(pi_decisions)

adata.uns["cell_type_final_v1_notes"] = {
    "leiden_resolution_used": 1.0,
    "method_basis":           "LLM verdict + marker dotplot",
    "rationale":              "LLM consensus high confidence on 14/18 clusters; PI revised cluster 7 (LLM said 'epithelial' but markers MUC5AC+TFF1+ → 'pit_cell')",
    "timestamp":              pd.Timestamp.now().isoformat(),
}
```

### Conventions

- **Key naming is PI's choice.** Use a meaningful prefix per analysis (`filter_v1`, `harmony_v1`, `gcpl_qc_2026_06_15`, etc.). Bumping to `_v2` follows the same versioning convention as stage h5ads.
- **No reserved framework namespace.** `adata.uns[key]` is yours; no `adata.uns["scrna_integration"][...]` to worry about.
- **The dict is plain Python.** Anything JSON-serialisable works (params, summary stats, cross-references to obsm/obs columns, free-text rationale).
- **Multi-version coexistence.** Old keys are never deleted on re-runs; instead a new `_v2` key records the rerun. This mirrors stage h5ad versioning so the same lineage rules apply inside one h5ad.

### When PI later wants to inspect what's been done

```python
[k for k in adata.uns if k.startswith("filter_") or k.startswith("harmony_")]
# → ['filter_v1', 'harmony_v1', 'harmony_v2']

adata.uns["harmony_v2"]
# {"method": "harmony", "batch_key": ["batch", "donor_id"], ...}
```

Plain dict access. No API to remember.

---

## obs Schema — Conventions, Not Enforcement

The framework documents recommended obs columns (the "three layers" below) but does **not** ship a runtime validator. PI inspects obs after `read_with_manifest` and judges whether the schema is OK for their analysis.

### Layer 1: Required core (manifest enforces, IO will fail loudly if missing)

`project_id`, `source_dataset`, `sample_id`, `donor_id`, `batch`, `cell_id` (auto-generated), `disease_system`.

`disease_system` is the broad research domain a project addresses (`"gastric"`, `"synovium"`, `"ovary"`, `"vaginal"`, ...). It is required because cross-disease integration analyses across PI's project portfolio (gastric / RA / PCOS / VVC) routinely group cells by this dimension; making it Layer 1 ensures it is always present and queryable. Per-project values are written via the manifest's top-level `disease_system` field.

### Layer 2: CellxGene-aligned (strong-warn + LLM best-effort fix)

`disease`, `disease_ontology_term_id`, `tissue`, `tissue_ontology_term_id`, `assay`, `sex`, `development_stage`.

All seven fields trigger a **strong warning** when missing or malformed at IO time. The framework attempts a best-effort fix using LLM-assisted disambiguation (looking at obs head + manifest + clinical-metadata head) and proposes a value; PI confirms; the confirmed fix is written back into `manifest.yaml` for reproducibility. Fields the LLM judges genuinely unfixable (e.g. a dataset that genuinely has no age information available) remain NaN, and PI sees the warning but is not blocked.

This applies uniformly to all seven fields rather than only the most-used ones, because each Layer 2 field is heavily used in at least one downstream stage:

- `disease` / `disease_ontology_term_id` — primary disease grouping in cross-condition DEG, abundance analysis, ontology-aware aggregation
- `tissue` / `tissue_ontology_term_id` — tissue-aware marker selection, 06 annotation, cross-tissue comparison
- `assay` — technical batch source for Harmony / scVI / integration QC
- `sex` — common confounder in cross-condition analyses across all disease systems
- `development_stage` — needed when comparing pediatric / adult / aged cohorts; relevant in PCOS (premenopausal vs postmenopausal) and gastric (age-correlated SPEM frequency)

### Layer 3: Project-defined (free)

Anything else: `disease_grade`, `treatment_arm`, `H_pylori_status`, etc. Pass through.

---

## Manifest Format

`data/{source_dataset}/manifest.yaml` is the only non-trivial schema in the framework. The manifest is the input to `read_with_manifest`. Whitelisted into git despite `data/` being gitignored.

### Write-once vs tunable: what belongs in the manifest

The manifest holds **write-once dataset facts** — how this dataset's obs columns are named, how its clinical table joins, which author-annotation columns exist, whether the author already removed doublets, species / ontology constants. Set once when a dataset is first ingested, then essentially frozen.

**Tunable analysis parameters** (QC thresholds, HVG count, clustering resolution, which embedding to promote) do **not** go in the manifest — they live in each stage notebook's `# === PARAMS ===` cell, where they are meant to be adjusted across re-runs. See ADR-0006 for why this split exists despite the 项目构思's general aversion to YAML.

Decision rule for any new config: *is it a fact about the dataset that's set once, or a knob that gets tuned?* Facts → manifest; knobs → PARAMS.

### Minimal manifest (the common case)

A clean dataset with symbol gene IDs and no clinical join needs only the six required fields — about eight lines:

```yaml
species: "human"
input:
  format: "10x_mtx"
  path: "filtered_feature_bc_matrix"
source_dataset: "GSE249874"
project_id: "gcpl_gastric_2026"
disease_system: "gastric"
original_annotations: []   # mandatory section; [] when the dataset ships no author labels
```

The seven optional blocks below (`obs_mapping`, `value_mapping`, `clinical_metadata`, `ontology`, `project_specific`, non-empty `original_annotations`, `qc_overrides`) appear **only when this dataset actually has the heterogeneity they describe**. Most datasets need two or three of them, not all seven.

### Full-featured manifest (reference superset — not a per-dataset requirement)

The example below shows every block at once for reference. No single dataset is expected to use all of them; copy only the blocks your dataset needs on top of the minimal six fields above.

```yaml
species: "human"               # mandatory; currently only "human" is accepted
                               # other values fail loudly (framework does not support cross-species)

input:
  format: "10x_mtx"            # 10x_mtx | h5ad | h5 | rds
  path: "filtered_feature_bc_matrix"
  raw_path: "raw_feature_bc_matrix"   # OPTIONAL; if present, written to adata.uns["raw_matrix_path"]
                                      # for 02 SoupX. 01 does NOT load it (avoids 2x memory).
                                      # Absent → SoupX skipped at 02.
  gene_id_format: "symbol"     # OPTIONAL; symbol | ensembl | auto (default auto)
                               # Hint to read_with_manifest about var.index format.
                               # In auto mode, framework detects format from var.index pattern.
  rds_object_class: "Seurat"   # only when format=rds
  rds_assay: "RNA"
  rds_slot: "counts"

source_dataset: "GSE249874"
project_id: "gcpl_gastric_2026"
disease_system: "gastric"      # mandatory; broad research domain (gastric | synovium | ovary | vaginal | ...)
                                # written into Layer 1 obs.disease_system for cross-project grouping

obs_mapping:
  sample_id: "Sample"
  donor_id:  "Patient"
  disease:   "Group"

value_mapping:
  disease:
    "CAG_mild":   "CAG"
    "CAG_severe": "CAG"
    "Healthy":    "normal"

ontology:
  disease_ontology_term_id: "MONDO:0005047"
  tissue:                   "gastric mucosa"
  assay:                    "10x 3' v3"

project_specific:
  # Layer 3 fields are entirely PI's choice; whatever subgrouping each project needs.
  # Below shows the common pattern (derive a categorical from an existing obs column)
  # and a non-exhaustive reference list of fields PI's projects have used to date.
  # Add / remove freely per project — these are examples, not framework requirements.
  disease_grade:
    source_column: "Group"
    rules:
      "CAG_mild":   "mild"
      "CAG_severe": "severe"
  # Other Layer 3 fields encountered across projects (illustrative, not enforced):
  #   timepoint            # baseline / week_4 / week_12 — any cohort subgroup, time or otherwise
  #   treatment_arm        # control / arm_A / arm_B
  #   H_pylori_status      # gastric projects: positive / negative / N/A
  #   OLGA_stage           # gastric pathology stage: 0 / I / II / III / IV
  #   OLGIM_stage          # gastric IM pathology stage
  #   DAS28                # RA disease activity score (continuous)
  #   remission_status     # RA / PCOS treatment-response
  #   PCOS_phenotype       # Rotterdam phenotype A / B / C / D
  #   BMI                  # continuous
  #   smoking_history      # never / former / current
  #   family_history       # disease-specific yes/no
  #   medication           # current medication categories
  #   surgery_type         # if relevant
  # New projects extend this list freely; the framework does not enforce a fixed schema here.

clinical_metadata:
  - file: "Supplementary_table_S1.xlsx"
    sheet: "Sheet1"
    skip_rows: 2
    join_on: {manifest_field: "sample_id", table_column: "Patient ID"}
    column_mapping:
      "Age": "age_years"
      "Sex": "sex"
      "H. pylori status": "H_pylori_status"
    value_mapping:
      sex:             {"M": "male", "F": "female"}
      H_pylori_status: {"+": "positive", "-": "negative", "N/A": null}
    on_missing:  "warn"           # warn | strict | silent
    on_conflict: "metadata_wins"  # metadata_wins | obs_wins | error

original_annotations:        # mandatory section, [] if dataset has no author annotations
  - column: "cell_type"
    role: "primary"
    granularity: "broad"
  - column: "celltype_subset"
    role: "subset"
    granularity: "fine"

preprocessing_done: ["basic_filter", "doublet_removal"]   # author-applied steps; QC skips matching ones

qc_overrides:
  doublet_removal:
    skip: true
    reason: "n_cells < 500, scrublet false positive rate too high"   # reason MANDATORY when skip:true
```

`read_with_manifest` validates the manifest schema (`species` present and equal to `"human"`; mandatory sections present; `qc_overrides[*].reason` provided when `skip:true`; `original_annotations` section present even if `[]`). Other validation is the caller's responsibility.

---

## R Bridge: subprocess Rscript for all R tools

R interoperability is unavoidable, but the mechanism is unified: **all R tools use subprocess `Rscript` + temp files**. See ADR-0007 (revised 2026-06-09). There is no rpy2 dependency in the pipeline. Three cases:

1. **Pure Python, no R** — InferCNV (`infercnvpy`: `cnv.tl.infercnv` / `cnv.tl.cnv_score`) and CytoTRACE (`cellrank` `CytoTRACEKernel`). These have mature Python packages validated in `student-code`; do not reach for R.
2. **subprocess `Rscript` + temp files** — SoupX, Monocle3, UCell, DESeq2, hdWGCNA, and all other R tools. Write `.mtx`/`.csv`, call `subprocess.run(["Rscript", "--vanilla", ...])`, read results back. Process isolation is robust; each `.R` script is independently debuggable; this is the pattern `student-code` actually uses throughout.
3. **rpy2 `%%R` in-notebook** — **Removed (2026-06-09)**. Formerly used for SoupX at 02. rpy2 + anndata2ri was infeasible in conda R 4.4.3 (`R_getVar` symbol missing); PI decision B moved SoupX to subprocess Rscript, unifying all R tools under one bridge mechanism.

The only framework-side R concession is that `read_with_manifest` handles `format: "rds"` internally (otherwise every IO call site needs rpy2 ceremony just to read a Seurat object). All other R access follows the subprocess pattern above, written directly in notebooks.

### SoupX subprocess pattern (02 ambient correction)

See `scripts/soupx_run.R` for the standalone R script. The notebook pattern:

```python
# Cell A: export filtered counts per sample (10x mtx format), extract original barcodes
# Cell B: call the independent R script
import subprocess
subprocess.run(["Rscript", "--vanilla", "scripts/soupx_run.R",
                work_dir, filtered_export_dir, raw_mtx_dir, sample_id], check=True)
# Cell C: read corrected counts back, update adata.X
```

The `.R` script (`SoupChannel` → `autoEstCont` → `adjustCounts` → write corrected `.mtx`) is a standalone file in `scripts/`, independently runnable and debuggable outside the notebook.

### Concrete example — SoupX (02 ambient correction, subprocess Rscript)

```python
# Cell A: R 环境守卫——检查 Rscript 可执行 + SoupX R 包可加载
import shutil, subprocess
_rscript = shutil.which(RSCRIPT_BIN)
_check = subprocess.run([_rscript, "--vanilla", "-e",
    'suppressPackageStartupMessages(library(SoupX)); cat("OK")'],
    capture_output=True, text=True, timeout=30)
r_available = _check.returncode == 0 and "OK" in _check.stdout

# Cell B: 对每个有 raw matrix 的样本导出过滤后计数，subprocess 调 Rscript
#   1. 提取原始 10x barcode（从 cell_id 前缀中剥离）
#   2. 导出过滤后矩阵（matrix.mtx + barcodes.tsv + features.tsv）到临时目录
#   3. subprocess.run([RSCRIPT_BIN, "--vanilla", "scripts/soupx_run.R", work, filt, raw, sid])
#   4. 读回 corrected_counts.mtx，写回 adata.X

# Cell C: 读回并写回
import scipy.io
corrected = scipy.io.mmread("results/_soupx_tmp/{sample}/corrected_counts.mtx")
adata[cell_ids].X = sp.csr_matrix(corrected.T)
adata.obs.loc[cell_ids, "ambient_correction_applied"] = True
```

### Concrete example — Monocle3 (10 pseudotime, subprocess Rscript)

Heavy R tools use subprocess + temp files (ADR-0007), the pattern validated in `student-code/workflow_for_pseudotime/4.3_*.py`:

```python
# Cell A: export from Python — sparse counts as .mtx (genes × cells), meta as .csv
import scipy.io, scipy.sparse as sp
work = "results/_monocle3_tmp"
scipy.io.mmwrite(f"{work}/counts.mtx", sp.csr_matrix(adata.X).T)   # R wants genes in rows
adata.obs.to_csv(f"{work}/cell_meta.csv")
adata.var.to_csv(f"{work}/gene_anno.csv")
# also export the chosen embedding (X_scVI / X_scANVI / X_pca) so Monocle3 reuses it

# Cell B: call an independent R script (lives at scripts/monocle3_trajectory.R)
import subprocess
subprocess.run(["Rscript", "--vanilla", "scripts/monocle3_trajectory.R", work], check=True)

# Cell C: read results back, write into adata + display the R-produced figures
import pandas as pd
pt = pd.read_csv(f"{work}/pseudotime.csv", index_col=0)
adata.obs["pseudotime_monocle3_v1"] = pt["pseudotime"].reindex(adata.obs_names)
# from IPython.display import Image; Image(f"{work}/trajectory_umap.png")
```

The `.R` script (`new_cell_data_set` → `preprocess_cds` → `reduce_dimension` → `cluster_cells` → `learn_graph` → `order_cells` → write `pseudotime.csv` + save plots) is a standalone file, independently runnable and debuggable outside the notebook.

### Concrete example — pseudobulk DESeq2 (08 pseudobulk DEG, subprocess Rscript)

```python
# Cell A: pseudobulk in Python (sample-level aggregation)
import decoupler as dc
pdata = dc.get_pseudobulk(adata, sample_col="sample_id", groups_col="cell_type_final_v1")

# Cell B: export counts + metadata, call DESeq2 script
work = "results/_deseq2_tmp"
pd.DataFrame(pdata.X.T, index=pdata.var_names, columns=pdata.obs_names).to_csv(f"{work}/counts.csv")
pdata.obs.to_csv(f"{work}/metadata.csv")
subprocess.run(["Rscript", "--vanilla", "scripts/deseq2_contrast.R",
                work, "disease", "CAG", "normal"], check=True)

# Cell C: read DEG table back
results = pd.read_csv(f"{work}/deg_CAG_vs_normal.csv", index_col=0)
adata.uns["deg_pseudobulk_CAG_vs_normal_v1"] = results.to_dict()
results.to_csv("results/deg_pseudobulk_CAG_vs_normal_v1.csv")
```

The `.R` script wraps `DESeqDataSetFromMatrix` → `DESeq` → `results(contrast=...)` → write CSV. Contrast levels are passed as `Rscript` args so the same script serves every contrast.

### Conventions for R cells

**subprocess (all R tools):**
- **The `.R` script lives in `scripts/`**, takes a working-directory path (and any tool-specific args) on the command line, and is independently runnable — debug it with `Rscript --vanilla scripts/foo.R results/_tmp` outside Jupyter.
- **Always pass `check=True`** to `subprocess.run` so a failing R script raises in Python instead of silently continuing.
- **Read R-produced figures back into the notebook** (`IPython.display.Image`) so results stay visible in the notebook even though computation happened out of process.
- **Temp dirs go under `results/_<tool>_tmp/`** (gitignored) and are cleaned at the end of the cell when large.
- **RSCRIPT_BIN convention**: each notebook PARAMS cell declares the Rscript path (e.g. `RSCRIPT_BIN = "~/miniforge3/envs/scrna-integration-r/bin/Rscript"`). Same convention across all R-using stages.

### rpy2: removed from the pipeline

rpy2 was formerly used only for SoupX at 02. The 2026-06-09 ADR-0007 revision reclassified SoupX to subprocess Rscript after rpy2 + anndata2ri proved infeasible in conda R 4.4.3 (`R_getVar` symbol missing). The rpy2 dependency is no longer needed for any stage. `read_with_manifest` retains rpy2 internally for `format: "rds"` (Seurat object reading), but this is the only remaining rpy2 call site in the codebase.

### Environment management

Two conda environment files at repo root: `environment.yml` (Python) and `environment-r.yml` (R). The R env is structured by stage so users skip optional packages (e.g. someone only doing QC skips Monocle3).

**环境隔离硬规定（agent 与 PI 都必须遵守）**：所有包安装一律在**专用命名 conda 环境**内进行，**绝不允许动 base / 主环境**。

- Python 环境固定命名 `scrna-integration`，R 环境固定命名 `scrna-integration-r`（写入 `environment.yml` / `environment-r.yml` 的 `name:` 字段）。
- agent **禁止** `pip install` / `conda install` 到 base 或任意非专用环境；禁止 `conda install -n base`；禁止无隔离的全局 `pip install`。
- 任何安装前必须先 `conda activate scrna-integration`（或 `-r`）；环境不存在则先 `conda env create -f environment.yml`（或 `-r.yml`）。
- CI 不实际装重型科学栈——`pyproject.toml` 与两个 `environment.yml` 是**声明**，真实安装由 PI/coder 在本地专用环境完成。验收以"依赖声明完整、smoke test（纯 import 框架包）通过、环境文件可被 conda 解析"为准，不要求 CI 跑完整安装。
- reviewer 审 env 相关 PR 时红线检查：diff 出现写向 base/全局/非专用环境的安装命令或文档指引 → flag。

### 跨平台一致性：精确 pin 源 spec + env_parity 诊断脚本 + 人工对齐（ADR-0010）

项目同时运行在两台机器：**开发机 Mac（Apple Silicon，`osx-arm64`）** 与 **服务器 Linux（x86-64，`linux-64`，无 GPU）**。要求两机的 conda 环境、包版本、代码函数行为**尽可能完全一致**；无法兼容的极少数包做**最小限度**显式登记切换。

**版本一致 — 精确 pin + 诊断 + 人工对齐**：

- **源 spec**（`environment.yml` / `environment-r.yml`）：人类可读，pin 到 Mac 已验证版本（`==`，**不带 build string** 以保持跨平台可移植）。这是"期望版本基准"的真相源。`pyproject.toml` 保持 `>=` 抽象依赖（库语义，不参与锁定）。
- **两机各自 `conda env create -f environment.yml`（或 mamba）直接安装**。精确 pin 已尽量保证一致。
- **诊断脚本 `scripts/env_parity.py`**（纯 Python 标准库 + subprocess 调 conda）：
  - `snapshot`：感知当前机器（`platform_tag()` 识别 linux-64/osx-arm64）、导出两个 conda 环境完整包清单（`conda list --json`）、写入 `docs/env-snapshots/{platform_tag}.json`（进 git）。
  - `compare`：对比两份快照（默认 linux-64 vs osx-arm64），输出版本不一致、仅 A 有、仅 B 有三类差异表。完全一致 exit 0，有差异 exit 1（便于 CI/hook 感知但不阻断）。
  - **硬约束：脚本只诊断报告差异，绝不自动修改环境**（对齐决策归人）。
- **人工对齐流程**：`compare` 出差异后，由人在落后那台执行 `conda install pkg=version` 或调 environment.yml 后重建。对齐后重新 `snapshot` 更新快照，提交 git 供下次对比。

**对齐异常登记** — `docs/cross-platform-exceptions.md`：

osx-arm64 的 bioconda 覆盖弱于 linux-64，少数 R/生信包可能缺 native build。这是允许"最小切换"的唯一合法场景：

- 优先用 conda-forge/bioconda 等价或 noarch build 对齐。
- 实在无法对齐的，**显式登记**：包名、缺失平台、根因、回退方案（Mac 走 Rosetta `osx-64` 子环境 / R 内 `BiocManager::install` / 仅某平台启用）、功能影响。
- 这是**异常清单非常态**，reviewer 每项质询"是否真无法对齐"，目标这张表尽可能短乃至为空。
- `monocle3` / `hdWGCNA` 本无 conda 包（走 R 内 `BiocManager`/`remotes`，两平台同一路径，天然一致），不算异常。

**代码一致 — OS 检测单点收口** `src/scrna_integration/platform.py`：

- 全项目**唯一**的 OS/路径差异收口点。`os.uname()` / `sys.platform` / `platform.system()` 判断**只允许**出现在此模块。
- `rscript_bin()`：从 `CONDA_PREFIX` 上溯派生 `{envs_dir}/scrna-integration-r/bin/Rscript`，跨平台同一逻辑。**所有 notebook 删除 Mac 硬编码 fallback 与写死 `"Rscript"`，统一改调 `platform.rscript_bin()`**。
- `platform_tag()`：返回标准化平台标识（`linux-64` / `osx-arm64` / `osx-64`），供快照文件命名与跨机器对比。
- 判据：notebook 与 `src` 其他文件**不得出现** OS 判断或平台绝对路径。reviewer grep `sys.platform` / `os.uname` / `/Users/` / `/home/` 硬编码路径，命中本模块以外即 flag。

**环境同步硬约束**：

- **conda 环境目录永不进 Syncthing / git**（环境内是平台相关编译二进制，arm64 `.so`/`.dylib` 在 linux-64 无法运行，跨平台同步会损坏环境）。只有源 spec、`platform.py`、异常登记表、快照 JSON 随同步走。两机各自 `conda env create` 安装。

**验收**：全流程需在 **Linux 与 Mac 两台机器各跑一遍端到端**，记录版本一致性 + 异常表，作为 ADR-0010 落地验收。

---

## Reference Data

Two read-mostly resources live under `references/` (whitelisted into git despite the rest of `references/` being gitignored):

### `references/markers/{tissue}_{purpose}.csv`

Per-tissue marker library. PI/students append rows over time. Schema:

```csv
tissue,cell_type,marker,role,reference,notes
gastric mucosa,SPEM,TFF2,canonical,Goldenring 2017,核心研究方向
gastric mucosa,SPEM,MUC6,canonical,Goldenring 2017,
gastric mucosa,pit_cell,MUC5AC,canonical,Nowicki-Osuch 2023,
```

Loaded via the framework's `load_markers` helper (see "The Three Functions" above for the full signature):

```python
from scrna_integration import load_markers

markers = load_markers("references/markers/gastric_epithelial.csv")
# → {"SPEM": ["TFF2", "MUC6", ...], "pit_cell": ["MUC5AC", ...]}
# Default role filter excludes 'negative'; pass roles=("negative",) or roles=None
# for other access patterns.
```

### `references/disease_ontology/{system}.yaml`

Per-system disease hierarchy. Loaded directly:

```python
import yaml
ontology = yaml.safe_load(open("references/disease_ontology/gastric.yaml"))
# PI walks the dict to compute ancestors / descendants in the analysis script
```

Schema (one node per entry, `parent` is null at roots):

```yaml
ontology: "gastric"
nodes:
  - {id: "CAG", label: "Chronic atrophic gastritis", parent: "chronic_gastritis", mondo: "MONDO:0005047"}
  - {id: "IM",  label: "Intestinal metaplasia",      parent: "CAG"}
  - {id: "IM_incomplete", label: "Incomplete IM",    parent: "IM"}
```

A minimal gastric ontology (5–10 nodes covering GCPL 5 datasets) ships with PR-1.

### What references/ does NOT contain (deliberately)

The following long-term-asset categories were considered and rejected:

| Considered | Rejected because |
|---|---|
| `references/cell_type_ontology/{system}.yaml` (cell type hierarchy) | LLM in 06 cross-method comparison already handles `"T cell"` vs `"CD8+ Tem"` semantic alignment without an explicit ontology. Maintaining a parallel cell type hierarchy adds cost without observed benefit. |
| `references/genesets/` (gene set library) | GSEApy / decoupler ship MSigDB / KEGG / Reactome built-in. Project-specific gene sets (e.g. SPEM transition signature) are best stored where they're used — in the relevant analysis notebook or per-project data dir — not promoted to a framework-wide library. |
| `references/embeddings/` (cellxgene_census ref metadata) | 04 notebooks call `cellxgene_census` API directly with hard-coded metadata. Locally archiving model versions adds maintenance without solving a real problem. |
| `references/prompts/` (LLM verdict / annotation prompts) | Per-stage LLM prompts live in the stage notebook scripts directly. Centralising them creates a sync burden between notebook and prompt file with no clear advantage. |

The decision rule for adding a third asset category to `references/` mirrors ADR-0001 / 0003 / 0004: a candidate must demonstrate observed real-world benefit (not anticipated convenience) and a new ADR before being added.

---

## QC Heterogeneity

QC heterogeneity is a first-class reality in multi-source scRNA-seq integration. Some source datasets ship pre-filtered (Nowicki-Osuch); some are organoids with inherently different MT% baselines (Kim, Yue); some are raw cellranger outputs (Nancang). A single set of QC thresholds applied blindly across all datasets is a scientific error.

The framework addresses this through **per-dataset QC notebooks** (ADR-0011). Each dataset's notebook under `notebooks/01_per_dataset/` independently applies the full QC pipeline with dataset-appropriate parameters. At merge time (`02_merged.ipynb`), per-dataset `adata.uns["qc_report_v1"]` records are aggregated into a cross-dataset summary, and `uns["merge_report_v1"]` records whether QC strategies were heterogeneous.

### Three QC strategies per dataset

Each per-dataset notebook declares its QC strategy in `adata.uns["qc_report_v1"]["strategy"]`:

1. **adaptive** — MAD-based thresholds (`median +/- N_MAD * MAD`) computed from the dataset's own distribution. The default; appropriate for most datasets.
2. **fixed** — explicit hard thresholds (e.g. `min_genes=200, max_genes=6000, max_pct_mt=20`). Used when PI has domain knowledge about a specific dataset's expected range.
3. **skip** — a QC step is skipped because the dataset author already performed it. Declared in manifest via `preprocessing_done` and `qc_overrides`. The notebook still column-aligns affected obs columns (fill NaN for skipped operations).

The MAD multiplier (`N_MAD`) is the real QC tuning knob — higher values admit more cells, lower values are stricter. The per-dataset notebook uses `N_MAD = 3` as default for `n_genes` and `total_counts`, and `N_MAD = 2.5` for `pct_counts_mt`.

### Mark-not-filter principle

All per-dataset QC follows the **mark-not-filter** principle: problematic cells are marked with boolean obs columns (`predicted_doublet`, `high_mt`, `high_ribo`, `stress_cell`, etc.) but not automatically removed. PI inspects the distributions (violins, scatter plots, per-dataset summaries) in the notebook and decides which cells to retain. This is a deliberate departure from the old approach of applying hard thresholds in batch — it gives PI the final say on filtering decisions that affect all downstream analyses.

### QC dimensions per dataset

Each per-dataset notebook runs the following QC dimensions:

1. **MAD-based adaptive filtering** — `n_genes`, `total_counts`, `pct_counts_mt` thresholds computed per dataset
2. **Doublet detection** — scrublet (skipped when manifest declares `preprocessing_done: ["doublet_removal"]`)
3. **Ambient RNA correction** — SoupX via subprocess Rscript (skipped when `input.raw_path` is absent)
4. **Cell cycle scoring** — `sc.tl.score_genes_cell_cycle` computes S/G2M scores and phase classification
5. **Gene complexity** — `n_genes_by_counts`, `log1p_total_counts` for richness-vs-depth diagnostics
6. **Special gene markers** — hemoglobin (HBA/HBB, erythrocyte contamination), stress response (JUN/FOS, dissociation artefact), ribosomal protein fraction

### Column alignment at merge

For downstream stages to work, QC-derived obs columns must be present and aligned across all integrated cells:

- `obs.n_genes`, `obs.total_counts`, `obs.pct_counts_mt`, `obs.pct_counts_ribo` — computed in per-dataset notebooks, preserved through merge.
- `obs.predicted_doublet`, `obs.doublet_score` — present as columns; **NaN** on cells where doublet detection was skipped (encodes "not applicable", not "missing data").
- `obs.ambient_correction_applied` (bool) — true on cells where SoupX ran; false elsewhere.
- `obs.phase` (S/G1/G2M) — from cell cycle scoring.
- `obs.high_mt`, `obs.high_ribo`, `obs.stress_cell` — mark-not-filter flags (PI reviews before actual removal).

### Cross-dataset QC diagnostics (Stage 02)

`02_merged.ipynb` provides cross-dataset QC diagnostics:

- **Violin plots** of `n_genes`, `total_counts`, `pct_counts_mt` per `source_dataset` — PI visually assesses whether distributions are broadly comparable
- **Gene intersection analysis** — bar charts of gene counts per dataset, union vs intersection size (join="inner" ensures all downstream stages operate on shared genes)
- **QC strategy summary table** — per-dataset `qc_report_v1` records aggregated into a single table showing strategy type, cells removed, and removal percentages
- **QC heterogeneity flag** — `uns["merge_report_v1"]["qc_heterogeneous"] = True` when datasets used different QC strategies; downstream stages can read this flag to include appropriate disclaimers

The key principle: cross-dataset QC diagnostics aim to **understand and document differences**, not to force alignment. If organoid datasets (Kim, Yue) have inherently higher MT% baselines than tissue biopsies (Nancang, Nowicki), forcing alignment would be a scientific error.

### What the framework does NOT enforce

PI explicitly accepts that QC heterogeneity propagates into 04-7 results:

- **Doublet alignment is NOT mandatory across datasets**. A dataset that ships already-doublet-removed (Nowicki) skips scrublet at 01; a raw cellranger dataset (Nancang) runs scrublet. The two datasets enter 03+ with different doublet histories. PI accepts the resulting risk that Harmony/scVI batch correction at 04 may treat residual doublets in some datasets as false batch effect. Mitigation is a disclaimer in the merge report when `qc_heterogeneous=True`.
- **Ambient correction is physically gated by data availability.** SoupX requires `input.raw_path` in the manifest. Datasets that ship only filtered matrices cannot run SoupX. `obs.ambient_correction_applied` records the actual situation per dataset.
- **Pseudobulk DEG / pseudotime** carries the same disclaimer. Pseudobulk averaging dilutes single-cell-level QC heterogeneity, but cannot eliminate systematic bias when one cohort has been doublet-removed and another has not.

### Cross-method comparison reads QC context

When 06 cross-method comparison includes original author annotations (`cell_type_original_{source_dataset}_v1`), the LLM verdict prompt is supplied with each source dataset's QC skip record (aggregated from per-dataset `qc_report_v1` entries) so the LLM can interpret author labels in their proper QC context. This is implemented in `06_annotated.ipynb` LLM verdict prompt construction — not in framework code.

---

## Memory Discipline (Conventions)

These are conventions documented here and enforced by **code review** (the `code-reviewer` agent flags violations). The framework does **not** ship runtime checks. Each convention has a concrete idiom PI / students copy verbatim into stage notebooks.

### 1. Sparse `adata.X` (CSR)

```python
import scipy.sparse as sp

# At ingest (already done by read_with_manifest, but explicit reminder for ad-hoc reads)
if not sp.issparse(adata.X):
    adata.X = sp.csr_matrix(adata.X)

# After any operation that might densify (some scanpy functions do this silently)
assert sp.issparse(adata.X), "adata.X became dense — investigate which step did it"
```

### 2. `inplace=True` by default

```python
sc.pp.filter_cells(adata, min_genes=200)         # modifies adata in place
sc.pp.normalize_total(adata, target_sum=1e4)     # modifies adata in place

# Only copy when you genuinely need a separate object (rare):
adata_before_qc = adata.copy()                   # explicit; PI knows memory just doubled
```

`adata.copy()` is the explicit signal that memory is now doubled. PRs that copy without justification get flagged.

### 3. Free at stage boundaries

```python
# Last cell of every stage notebook:
adata.write_h5ad(output_path, compression="lzf")
print(f"Wrote {output_path}")

del adata
import gc
gc.collect()
```

Without this, the Jupyter kernel keeps the previous stage's adata alive when PI runs the next stage in the same kernel. Across 7 stages × multiple versions on a 25k-cell pilot, omitting this can OOM a 16GB Mac.

### 4. `compression="lzf"` on every write

```python
adata.write_h5ad(output_path, compression="lzf")
```

Always pass it. `legacy-GCPL` notebooks omit this and the result h5ads run ~30% larger than necessary.

### 5. `float32` dtype throughout

```python
# After normalize:
adata.X = adata.X.astype(np.float32)

# When scvi/scanvi outputs land in obsm:
adata.obsm["X_scvi"] = adata.obsm["X_scvi"].astype(np.float32)
adata.obsm["X_scanvi"] = adata.obsm["X_scanvi"].astype(np.float32)
```

scRNA-seq counts and embeddings carry far less than float64's 15 decimal digits of precision. float64 doubles memory for no useful precision.

### In-notebook self-check (one-line assertion)

Every stage notebook places a single-line assertion immediately before `adata.write_h5ad(...)`. This catches the highest-impact memory regression — `adata.X` going dense or losing float32 — at the moment it would otherwise leak into the next stage's checkpoint:

```python
# Memory discipline self-check (one cell before write)
import scipy.sparse as sp
assert sp.issparse(adata.X) and adata.X.dtype == np.float32, \
    f"adata.X invariants violated: sparse={sp.issparse(adata.X)}, dtype={adata.X.dtype}"
```

Why this rather than a fuller check:

- `adata.X` is the largest object in memory (cells × genes). Guarding it captures ~95% of memory regression risk. `obsm` latent matrices (PCA / Harmony / scVI / scANVI ≈ 4–5 entries × cells × ~30–50 dims) are an order of magnitude smaller — guarding them is not worth a `for` loop / extra cell.
- Single-line assertion runs in nanoseconds; no measurable overhead at write time.
- `obsm` dtype regressions are caught at PR review time via the Reviewer Cheat Sheet (ADR-0004) — code-review-side defence, not runtime-side.

### Reviewer enforcement

When `code-reviewer` agent reviews a PR with framework or notebook code, it flags:

- `adata.X` becoming dense (hard violation)
- `.copy()` calls without justification in surrounding code or PR description (medium issue)
- Stage notebooks that don't end with the `del adata; gc.collect()` idiom (medium issue)
- `adata.write_h5ad(...)` calls without `compression="lzf"` (low issue)
- `float64` dtype on `adata.X` post-normalize or on `obsm/X_*` matrices (low issue)

Severity levels follow `code-reviewer` agent's standard verdict scale (block / request-changes / approve). Memory violations rarely block — they are corrections requested before merge.

---

## 06: Multi-method Annotation in the Notebook

06 is the most distinctive scientific layer of the framework, but it is implemented as a **notebook** (`06_annotated.ipynb`), not a framework module. Annotation methods run as **independent cells in the notebook** — PI is free to choose any execution order; the framework does not anchor "correct order" between manual marker annotation and the automatic methods. (Order experimentation is encouraged: running automatic methods first lets PI cross-check intuition, while running manual marker first preserves PI's prior judgement free of automatic-method anchoring.)

**Running several methods in parallel is the point, not redundancy.** The most accurate per-cluster cell type comes from cross-comparing independent methods — agreement raises confidence, disagreement flags clusters needing PI attention. So the default is to run multiple methods together, not to pick one.

Methods fall into three tiers:

**Default co-run (run together every time):**

1. **Marker dotplot** — PI views dotplots of canonical/optional markers from `references/markers/{tissue}_{purpose}.csv` and manually labels each leiden cluster, writing `obs.cell_type_marker_v1`.
2. **mLLMCelltype consensus via provider-direct API keys** — writes `obs.cell_type_llm_v1` + `_uncertainty` + `_consensus_log`. Same code pattern as `legacy-GCPL/06_annotation.ipynb`. Uses per-provider API keys + base URLs (no OpenRouter middleware by default). Model selection in `models=[...]` list; API routing by model name prefix. See `.env.example` for key/URL template.
3. **Gene-set scoring** — AUCell / UCell / scanpy `score_genes` over canonical marker sets, producing continuous scores in `obs.score_*` columns. Cross-cluster summarisation (mean score + pct-positive per cluster, transition/mixed flagging) follows the pattern in `student-code/.../6.3_UCell_cluster_mean_score&pct_pos.py`. Used as supplementary evidence in cross-method comparison rather than as standalone cell-type labels.

**Default-active when a reference atlas is available (data-availability gated, not preference):**

4. **scANVI label transfer** — writes `obs.cell_type_scanvi_v1` + `_uncertainty` from a labelled reference atlas. The cell is active by default; it runs whenever a suitable reference atlas exists for the disease system. When none is available the cell is skipped (it is the only default method with a hard external prerequisite — a labelled atlas for this tissue). The notebook cell header states the prerequisite explicitly.

**Candidate (commented-out cell, enable when applicable):**

5. **CellTypist pretrained classifier** — writes `obs.cell_type_celltypist_v1` + `_uncertainty`. Left as a commented-out cell with a header noting "enable when a CellTypist model exists for this tissue." Not run by default because PI's gastric / synovial systems do not currently have well-matched pretrained models; enabling it is a one-line uncomment when a relevant model appears.

After the method cells run, the notebook constructs the cross-method comparison over whichever methods produced labels:

6. **Confusion matrices + Sankey diagrams** between method pairs (with Cohen's kappa) using scanpy / seaborn natively.
7. **LLM verdict per cluster** — the notebook constructs a prompt containing all available method labels (including original-author annotations from `cell_type_original_*` if present), top markers from `sc.tl.rank_genes_groups`, gene-set score profile, marker_db coverage, and `adata.uns["qc_skipped"]` per source dataset. OpenRouter is called directly. The verdict for each cluster is written to a per-cluster markdown file.
8. **PI reads every cluster's verdict and decides** — `obs.cell_type_final_v1` is set manually. With one integrated adata per project (per-project N is typically 15–30 leiden clusters), reading every verdict takes 30–60 minutes and is the standard workflow. There is no auto-pick by confidence threshold — every cluster's final label is PI's call.

**There is no `si.report.06_annotation` function.** The notebook is the unit of work. When new annotation methods appear, PI / students add a cell to `06_annotated.ipynb` and rename the obs column following `cell_type_{method}_v{N}` naming.

### `cell_type_final_v{N}` versioning

PI may revisit annotations weeks later (e.g. after seeing a 07 result that suggests refining cluster 7 from `SPEM` into `SPEM_complete` / `SPEM_incomplete`). Following the same convention as stage h5ad versioning, **the previous `cell_type_final_v{N}` column is never overwritten** — a new `cell_type_final_v2` column is added with the revised labels, and `adata.uns["cell_type_final_v2_notes"]` records the revision rationale. PI / agents reference whichever version is appropriate for the analysis at hand.

### 06c subset re-annotation reflows back

When 06c produces refined sub-cluster annotations (e.g. CD4/CD8/Treg/MAIT from a T-cell subset re-clustering), the refined labels are reflowed into the main 06 adata as a new column `cell_type_final_subset_v1`. Cells in the refined subset receive the fine label; cells outside the subset (epithelial, stromal, etc.) receive `NaN` in this column. The original `cell_type_final_v1` column remains untouched. This preserves the granularity hierarchy (broad in `_final_v1`, fine in `_final_subset_v1`) without polluting either, and downstream stages can choose the column appropriate to their question.

### Original author annotations

Manifest's `original_annotations` section causes IO to rename author labels to `cell_type_original_{source_dataset}_v1[_{role}]`. After multi-source integration these columns are sparse (NaN on cells from datasets without that annotation). The 06 notebook treats them as additional methods in the cross-method comparison — the LLM verdict prompt naturally references them when present and reads `adata.uns["qc_skipped"]` so the LLM can interpret author labels in their proper QC context (see "Cross-method comparison reads QC context" in the QC Heterogeneity section).

### Sweep recommendations (no longer bound to `sweep()`)

The 06 notebook ends with one final cell that constructs an LLM call to generate **sweep recommendations** based on the single-run results: which clusters had high uncertainty across methods, which methods showed systematic disagreement, which marker_db gaps the LLM observed. Output is appended to the 06 markdown report. PI reads recommendations and decides whether to re-run annotation with different parameters (e.g. different LLM model sets, different scANVI references). Cost is one LLM call per 06 run — acceptable. This is the operational form of the project's "teach AI to make initial scRNA-seq judgements" thesis.
*Note: the `sweep()` function was removed per ADR-0009; the term "sweep recommendations" is a retained workflow concept, not a framework function.*

### Per-cluster deep profile

A second notebook (`06b_per_cluster.ipynb`) loops over clusters and emits one markdown per cluster with UMAP highlight, top-marker dotplot, gene-set score violins, cross-disease abundance boxplot, and an LLM-written narrative paragraph. **No framework class, no panel registry — just a notebook with a `for cluster in adata.obs[label_col].cat.categories: ...` loop.**

---

## Notebooks (Directly Runnable, Not Templates)

The `notebooks/` directory contains **directly runnable notebooks**, not a template-instantiation system. PI edits the PARAMS cell at the top and runs all cells — no copy / rename ceremony for Stage 02 onwards. Stage 01 is the one exception: it ships four notebooks named by input data format (`01_template_*`), each runnable as-is against the `data/_subset/` fixtures and intended to be copied per dataset. The framework's standardisation ships as the actual notebook content, not as a generator.

```
notebooks/
├── 01_per_dataset/              # per-dataset independent QC (ADR-0011)
│   ├── 01_template_10x_mtx.ipynb        #   10x mtx (CellRanger filtered_feature_bc_matrix/) — adaptive QC + scrublet + SoupX + cell cycle
│   ├── 01_template_10x_h5.ipynb         #   10x h5 (filtered_feature_bc_matrix.h5) — read_10x_h5 + var_names_make_unique
│   ├── 01_template_h5ad.ipynb           #   pre-processed h5ad (CELLxGENE / HCA) — .raw extraction + obs_mapping + skip mode
│   └── 01_template_counts_matrix.ipynb  #   tsv_matrix (gzip tab-separated counts table) — manual table read + QC metric back-fill
├── 02_merged.ipynb              # anndata.concat(join="inner") + cross-dataset diagnostics + QC report aggregation
├── 03_normalized.ipynb          # normalize + log1p (+ optional Pearson residuals) + batch-aware HVG + HVG exclusion list + scalar-or-sweep
├── 04_embedded.ipynb            # PCA + Harmony + scVI + scANVI; elbow plot + N_NEIGHBORS sweep + HARMONY_THETA sweep + integration metrics
├── 05_clustered.ipynb           # multi-resolution Leiden + clustering stability (subsample ARI) + marker gene preview
├── 06_annotated.ipynb           # 5-method annotation + cross-method comparison + LLM verdict
├── 06b_per_cluster.ipynb        # per-cluster deep profile
├── 06c_subset.ipynb             # subset re-cluster (T cells / epithelial / etc.)
├── _deprecated/                 # superseded by per-dataset + merge design (ADR-0011)
│   ├── 01_loaded.ipynb          #   (deprecated) old monolithic multi-dataset load
│   └── 02_qcd.ipynb             #   (deprecated) old monolithic post-merge QC
└── 07_downstream/               # downstream modules; each consumes 06 (or 06c) output
    ├── deg.ipynb                    # PR-3 — sc.tl.rank_genes_groups
    ├── pseudobulk_deg.ipynb         # PR-3 — DESeq2 (subprocess Rscript)
    ├── cnv.ipynb                    # PR-3 — InferCNV (pure Python infercnvpy)
    ├── pathway.ipynb                # PR-4+
    ├── pseudotime.ipynb             # PR-4+
    ├── grn.ipynb                    # PR-4+
    ├── cell_communication.ipynb     # PR-4+
    ├── abundance.ipynb              # PR-4+
    └── gene_modules.ipynb           # PR-4+
```

### Common notebook structure

Every notebook follows the same structural template (not a code template — a structural pattern):

1. **PARAMS cell** at the top (`# === PARAMS ===` markdown header + Python cell with assignments).
2. **Imports cell** (scanpy / framework imports / subprocess + shutil if R tools used).
3. **Load upstream h5ad** (or call `read_with_manifest` for 01).
4. **Stage logic cells** — scanpy native APIs.
5. **Run-metadata cells** — `adata.uns[...] = {...}` writes documenting what was done.
6. **Report cells** — `sc.pl.*` + `plt.savefig` + per-stage markdown summary.
7. **Final cell** — `adata.write_h5ad(out_path, compression="lzf")` + `del adata; gc.collect()`.

### Per-notebook cell sequence (specification for PR-3 implementation)

The cell sequences below are the **specification** PR-3 coder agents implement against. Each line is one notebook cell; markdown cells are prefixed `[md]`, code cells are `[code]`.

#### `01_per_dataset/{dataset}.ipynb` (per-dataset independent QC)

One notebook per source dataset. Each follows the same structural template with dataset-specific thresholds in the PARAMS cell. The notebook is the QC report — all diagnostic plots are output directly inline; there is no separate report layer.

```
[md]   # 01: Per-dataset ingest + QC — {dataset_name}
[code] # === PARAMS ===
       MANIFEST_PATH = "data/{source_dataset}/manifest.yaml"
       OUTPUT_PATH   = "results/01_{dataset}_v1.h5ad"
       N_MAD_n_genes     = 3      # MAD multiplier for n_genes threshold
       N_MAD_total_counts = 3     # MAD multiplier for total_counts threshold
       N_MAD_pct_mt       = 2.5   # MAD multiplier for pct_mt threshold
       ENABLE_SCRUBLET    = True  # False when manifest declares doublet_removal already done
       ENABLE_SOUPX       = True  # False when manifest lacks raw_path
       RSCRIPT_BIN = None         # set to Rscript path or use platform.rscript_bin()
       SOUPX_RSCRIPT = "scripts/soupx_run.R"
       RANDOM_SEED = 42
[code] # imports (scanpy, scrublet, scipy, subprocess, platform)
[code] # call read_with_manifest, inspect returned adata
[md]   ## Raw data overview
[code] # print obs/var head, baseline metrics summary, obsm/uns keys
[md]   ## QC visualisation pre-filter
[code] # violin + scatter of n_genes / total_counts / pct_mt / pct_ribo
       # mark-not-filter: MAD thresholds computed, high_mt/high_ribo flags set
[md]   ## Doublet detection (conditional)
[code] # if manifest declares doublet_removal already done:
       #     obs.doublet_score = NaN, predicted_doublet = False, log skip reason
       # else:
       #     scrublet per sample (10x default params), store score + boolean
       #     plot doublet score histogram
[md]   ## Ambient RNA correction via SoupX (conditional)
[code] # guard 1: check adata.uns["raw_matrix_path"]
       # guard 2: check Rscript available + SoupX R package loadable
       # for each sample with raw matrix: subprocess Rscript → read corrected counts
       # if skipped: obs.ambient_correction_applied = False, log reason
[md]   ## Cell cycle scoring
[code] # sc.tl.score_genes_cell_cycle → obs.S_score, G2M_score, phase
       # violin of S/G2M scores
[md]   ## Gene complexity
[code] # sc.pp.calculate_qc_metrics with log1p=True → n_genes_by_counts, log1p_total_counts
       # scatter: n_genes vs total_counts coloured by pct_mt
[md]   ## Special gene markers
[code] # compute pct HBA/HBB (erythrocyte), JUN/FOS (stress), ribosomal fraction
       # mark-not-filter: stress_cell, high_ribo flags set
[md]   ## QC report & column alignment
[code] # obs column cleanup: ensure all expected columns exist (NaN-fill for skipped ops)
       # write adata.uns["qc_report_v1"] = {strategy, n_cells_in, n_cells_out, cells_removed,
       #   pct_removed, steps_skipped, skip_reasons, N_MAD_values, timestamp}
       # adata.uns["stage"] = "01_per_dataset"
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write h5ad, del adata, gc.collect
```

**Convention**: the four Stage-01 notebooks are named by **input data format**, not by dataset — a dataset belongs to the downstream project, whereas the format is the reusable skeleton the framework should supply. When adding a new dataset, copy the template matching its input format (`01_template_10x_mtx` for a CellRanger `filtered_feature_bc_matrix/` directory, `01_template_10x_h5` for a `.h5`, `01_template_h5ad` for an author-preprocessed `.h5ad`, `01_template_counts_matrix` for a gzip tab-separated counts table), rename the copy after the dataset, fill in the placeholders in the PARAMS cell (manifest path, `DATA_SOURCE_ID`, dataset-specific thresholds), and add the output path to `02_merged.ipynb`'s `PER_DATASET_PATHS` list. The templates ship pointing at the `data/_subset/` fixtures, so a fresh checkout runs as-is.

#### `02_merged.ipynb` (cross-dataset merge + diagnostics)

```
[md]   # 02: Multi-dataset merge + cross-dataset QC diagnostics
[code] # === PARAMS ===
       PER_DATASET_PATHS = [
           "results/01_<dataset_a>_v1.h5ad",
           "results/01_<dataset_b>_v1.h5ad",
       ]
       OUTPUT_PATH = "results/02_merged_v1.h5ad"
       JOIN_GENES = "inner"
       MIN_SHARED_GENES = 15000
       BATCH_KEY = "source_dataset"
       DOWNSAMPLE_TO_MIN = False
       OUTPUT_VERSION = 1
       RANDOM_SEED = 42
[code] # imports (scanpy, anndata, numpy, pandas, scipy, matplotlib)
[md]   ## Load per-dataset h5ads
[code] # for each path: sc.read_h5ad, collect adatas list + dataset_info table
       # display summary: dataset name, n_cells, n_genes, QC strategy, cells_removed
[md]   ## Gene intersection analysis
[code] # compute gene union/intersection across all datasets
       # bar chart: per-dataset gene counts + union vs intersection
       # assert n_shared >= MIN_SHARED_GENES
[md]   ## Merge (anndata.concat)
[code] # anndata.concat(adatas, join=JOIN_GENES, label="source_dataset", keys=src_keys, index_unique="-")
       # validate cell_id uniqueness
       # free per-dataset adatas; gc.collect()
[md]   ## Cross-dataset QC diagnostics
[code] # violin plots: n_genes / total_counts / pct_counts_mt per source_dataset
       # summary table: per source_dataset median genes, median UMI, median MT%, doublet_rate, phase distribution
[md]   ## QC heterogeneity record
[code] # aggregate per-dataset qc_report_v1 into uns["merge_report_v1"]
       # fields: n_datasets, join_genes, n_shared_genes, gene_intersection_pct,
       #          qc_heterogeneous (bool), qc_strategies (list), per_dataset_qc (dict)
[md]   ## Optional downsampling
[code] # if DOWNSAMPLE_TO_MIN: downsample each dataset to min cell count
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write h5ad, del adata, gc.collect
```

#### `03_normalized.ipynb`

```
[md]   # 03: normalize + log + HVG
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, n_top_genes, hvg_flavor
[code] # imports, load
[code] # adata.layers["counts"] = adata.X.copy()  (preserve raw counts before normalize)
[code] # sc.pp.normalize_total(adata, target_sum=1e4)
[code] # sc.pp.log1p(adata)
[code] # sc.pp.highly_variable_genes(adata, n_top_genes=..., flavor=...)
[code] # adata.X = adata.X.astype(np.float32)  (memory discipline)
[code] # sc.pl.highly_variable_genes(adata) + plt.savefig
[code] # adata.uns["normalize_v1"] = {target_sum, hvg_flavor, n_top_genes, timestamp}
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write, del, gc
```

#### `04_embedded.ipynb`

```
[md]   # 04: multi-method embedding (all peers) + integration sweep
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, n_pcs, embedding_methods (list), batch_key
[code] # imports + load
[md]   ## PCA (baseline)
[code] # sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)  → obsm["X_pca"]
[md]   ## Harmony
[code] # sce.pp.harmony_integrate(adata, key=batch_key)  → obsm["X_pca_harmony"]
[md]   ## scVI (train from scratch)
[code] # scvi.model.SCVI setup + train  → obsm["X_scVI"]
[md]   ## (Optional) cellxgene_census pretrained scVI
[code] # commented-out / conditional: download census model, mygene symbol↔ensembl align,
       #   prepare_query_anndata, load_query_data, get_latent_representation → obsm["X_scVI_census"]
       #   PREREQUISITE: census has a good pretrained model for this tissue
       #   (see legacy-GCPL/04_dimensionality_reduction.ipynb for the full pattern)
[md]   ## (Optional) scANVI label transfer embedding
[code] # conditional on a labelled reference: → obsm["X_scANVI"]
[md]   ## Visual comparison (PRIMARY decision input)
[code] # for each obsm: sc.pl.umap coloured by sample_id, batch, cell_type
       #   PI eyeballs batch integration vs biological over-mixing
[md]   ## Integration sweep — explicit for loop across embeddings with integration_metrics
[code] # for rep in [k for k in ["X_pca", "X_pca_harmony", ...] if k in adata.obsm]:
       #     adata_copy = adata.copy()
       #     sc.pp.neighbors(adata_copy, use_rep=rep)
       #     sc.tl.umap(adata_copy)
       #     m = integration_metrics(adata_copy)
       #     results.append({"use_rep": rep, **m})
       # sweep_df = pd.DataFrame(results); write report; display table
[code] # PI weighs UMAP visuals + metric table, decides which to mark promoted
[code] # adata.uns["harmony_v1"] / ["scvi_v1"] / etc. metadata writes
       # adata.uns["status"] = "experimental" (PI manually changes to "promoted" later)
[code] # cast all obsm latent matrices to float32
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write, del, gc
```

**Adding an embedding method**: write `adata.obsm["X_{method}"]` and add `"X_{method}"` to the sweep `use_rep` candidates — one cell, no framework change (see "04 sweep includes integration QC" above).
[code] # adata.uns["harmony_v1"] / ["scvi_v1"] / ["scanvi_v1"] metadata writes
       # adata.uns["status"] = "experimental" (PI manually changes to "promoted" later)
[code] # cast all obsm latent matrices to float32
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write, del, gc
```

#### `05_clustered.ipynb`

```
[md]   # 05: multi-resolution Leiden + sweep
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, resolutions (list), use_rep
[code] # imports + load
[md]   ## Multi-resolution Leiden
[code] # for res in resolutions: sc.tl.leiden(adata, resolution=res, key_added=f"leiden_res_{res}")
[md]   ## Clustering sweep — explicit for loop across resolutions with clustering_metrics
[code] # for res in resolutions:
       #     adata_copy = adata.copy()
       #     sc.tl.leiden(adata_copy, resolution=res, key_added=f"leiden_res_{res}")
       #     m = clustering_metrics(adata_copy)
       #     results.append({"resolution": res, **m})
       # sweep_df = pd.DataFrame(results); write report; display table
[md]   ## (Optional) alternative clustering method — see "Adding a clustering method" below
[code] # commented-out slot: any method that writes obs["{method}_clusters"] coexists
       #   with the leiden_res_* columns and is compared the same way
[md]   ## Visualisation
[code] # sc.pl.umap with each leiden_res_* column for PI to compare
[code] # uns metadata writes
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write, del, gc
```

**Adding a clustering method (extension pattern).** Multi-resolution Leiden is the default because it is scanpy-native and pairs with `sweep` to let PI pick the resolution from a scored table (consistent with SOUL's "judgement not outsourced"). But Leiden alone sometimes leaves the cluster count ambiguous, so the project must stay open to plugging in other methods (ACDC, future community-detection methods, etc.) as they appear — per the 项目构思's "学术追新" goal.

No framework code is needed for this. The extension mechanism is the same scanpy-native "multiple results coexist as parallel slots" pattern as 04 embeddings:

- Any clustering method writes its labels to its own obs column: `obs["acdc_clusters"]`, `obs["{method}_clusters"]`, etc. Multi-resolution Leiden already does this (`leiden_res_0.5`, `leiden_res_1.0`, ...).
- Methods that sweep a parameter (Leiden over resolution) go through `sweep()`. Methods that auto-pick (ACDC searches for an optimal partition itself) just write their single result column directly — no sweep needed.
- All resulting cluster columns are compared the same way: `sc.pl.umap` coloured by each, plus ARI / silhouette across columns. PI picks which column feeds 06.
- A new method is one added cell, not a framework change. ACDC specifically was tried in GCPL but never completed a run (prohibitively slow on PI's data), so it is **not** a default and `acdc_py` is **not** a PR-0a dependency; it is one example of the kind of method this slot accepts when PI chooses to install it.

#### `06_annotated.ipynb`

```
[md]   # 06: Multi-method annotation + cross-method comparison
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, MARKER_CSV, leiden_col
       # API keys loaded from .env (see .env.example for template)
[code] # imports + load + load_markers(MARKER_CSV)
[md]   ## Method 1: marker dotplot (PI manually labels)
[code] # sc.pl.dotplot(adata, var_names=markers, groupby=leiden_col)
[code] # PI inspects, fills marker_assignments dict, then:
       # adata.obs["cell_type_marker_v1"] = adata.obs[leiden_col].map(marker_assignments)
[md]   ## Method 2: mLLMCelltype consensus via provider-direct API keys
[code] # exact pattern from legacy-GCPL/06_annotation.ipynb
[md]   ## Method 3: gene-set scoring (AUCell / UCell / scanpy_score)
[code] # for cell_type, marker_list in markers.items(): sc.tl.score_genes(...)
       # then per-cluster mean score + pct_pos summary (student-code 6.3 pattern)
[md]   ## Method 4: scANVI label transfer (active when a labelled reference atlas exists)
[code] # PREREQUISITE: labelled reference atlas for this disease_system.
       # if REFERENCE_ATLAS_PATH set: run scANVI transfer → cell_type_scanvi_v1 + _uncertainty
       # else: skip (print why); cross-method comparison proceeds without it
[md]   ## Method 5 (candidate, commented out): CellTypist
[code] # ENABLE when a CellTypist model exists for this tissue:
       # # import celltypist
       # # adata.obs["cell_type_celltypist_v1"] = celltypist.annotate(adata, model="...").predicted_labels
[md]   ## Cross-method comparison: confusion matrices + Sankey
[code] # for each available method pair: confusion matrix + Cohen's kappa + sankey
[md]   ## LLM verdict per cluster
[code] # construct prompt with all method labels + top markers + scores + qc_skipped context
       # call LLM via provider API key (from .env), write per-cluster verdicts to results/figures/06_verdicts/
[md]   ## PI sets cell_type_final_v1 (manual review against verdicts)
[code] # pi_decisions = {leiden_id: cell_type, ...}
       # adata.obs["cell_type_final_v1"] = adata.obs[leiden_col].map(pi_decisions)
       # adata.uns["cell_type_final_v1_notes"] = {...}
[md]   ## Sweep recommendations (LLM-generated)
[code] # one OpenRouter call producing recommendations markdown appended to 06 report
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write, del, gc
```

#### `06b_per_cluster.ipynb`

```
[md]   # Per-cluster deep profile (one markdown per cluster)
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_DIR, label_col, geneset_panels (list of csv paths)
[code] # imports + load
[code] # for cluster in adata.obs[label_col].cat.categories:
       #     # subset, plot UMAP highlight, top markers dotplot, score violins,
       #     # cross-disease abundance boxplot, custom panels (list of plain functions),
       #     # LLM-written narrative paragraph (one OpenRouter call)
       #     # write cluster_{cluster}.md to OUTPUT_DIR
[code] # write index.md linking all cluster files with one-line summaries
```

#### `06c_subset.ipynb`

```
[md]   # 06c: subset re-cluster (e.g. T cells, SPEM-spectrum, epithelial)
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, subset_filter (str expression on obs), n_top_genes, resolutions
[code] # imports + load
[code] # adata_sub = adata[adata.obs.eval(subset_filter)].copy()
       # adata_sub.uns["subset_of"] = UPSTREAM_PATH
       # adata_sub.uns["subset_filter"] = subset_filter
[code] # re-run 03 (HVG re-selection on subset)
[code] # re-run 04 (re-embedding)
[code] # re-run 05 (re-clustering)
[code] # re-run 06 (re-annotation; usually finer cell types)
[code] # reflow refined labels back to main adata as cell_type_final_subset_v1 column
       # main_adata = sc.read_h5ad(UPSTREAM_PATH)
       # main_adata.obs["cell_type_final_subset_v1"] = ...  # fill subset cells, NaN elsewhere
       # main_adata.write_h5ad(UPSTREAM_PATH_with_subset_v1)  # bumped version of main h5ad
[code] # write subset h5ad + write updated main h5ad
```

#### `07/deg.ipynb` (PR-3 in scope)

```
[md]   # 07: per-cluster differential expression
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, label_col, method
[code] # imports + load
[code] # sc.tl.rank_genes_groups(adata, groupby=label_col, method="wilcoxon")
[code] # sc.pl.rank_genes_groups + sc.pl.rank_genes_groups_dotplot
[code] # export top-N DEG per cluster as CSV
[code] # adata.uns["rank_genes_groups_v1"] is set by scanpy; record additional metadata
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write, del, gc
```

#### `07/pseudobulk_deg.ipynb` (PR-3 in scope)

```
[md]   # 07: pseudobulk cross-condition DEG via DESeq2 (subprocess Rscript)
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, sample_col, contrast_col, contrast_levels
[code] # imports (decoupler for pseudobulk, subprocess for DESeq2 Rscript)
[code] # pseudobulk per (sample_id × cell_type_final_v1) using decoupler
[code] # export counts + metadata → subprocess.run([Rscript, "scripts/deseq2_contrast.R", ...])
[code] # collate results across cell types into one CSV / dict in adata.uns
[code] # plot: heatmap of top DEGs per cell type, volcano per contrast
[code] # write CSVs to results/, optionally adata.uns["pseudobulk_deg_v1"]
```

#### `07/cnv.ipynb` (PR-3 in scope)

```
[md]   # 07: CNV inference via infercnvpy (pure Python, no R)
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, reference_cells_col, reference_cats, gene_position_file
[code] # imports (import infercnvpy as cnv) + load
[code] # add chromosomal positions to adata.var (chromosome / start / end);
       #   sort by chromosome+position; drop chromosomes with <20 genes
       #   (see student-code 4.2_add_gene_positions_for_infercnv.py)
[code] # cnv.tl.infercnv(adata, reference_key=reference_cells_col, reference_cat=reference_cats,
       #                 window_size=100, step=10, lfc_clip=3.0)
[code] # cnv.tl.pca(adata); cnv.pp.neighbors(adata); cnv.tl.leiden(adata)
[code] # cnv.tl.cnv_score(adata, groupby="cnv_leiden")  → adata.obs["cnv_score"]
[code] # cnv.pl.chromosome_heatmap; identify malignant vs reference clusters
[code] # adata.uns["infercnv_v1"] metadata
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write, del, gc
```

### Conventions across all notebooks

- Every notebook starts with `# === PARAMS ===` cell — uppercase param names, explicit string paths, no magic resolution.
- Every notebook ends with the memory discipline self-check (one-line assertion on `adata.X`) followed by `adata.write_h5ad(out, compression="lzf")` + `del adata; gc.collect()`.
- Every notebook uses scanpy native APIs and the two framework functions (`read_with_manifest` / `load_markers`) plus scorers (`from scrna_integration.scorers import ...`) only where they fill a real gap.
- Run-metadata writes follow the "Run Metadata" section above (plain `adata.uns[key] = {...}` with versioned keys).
- **Cell-sequence specs are binding.** When a notebook drifts from its spec (cells removed / reordered / renamed) during PR review, that is a request-changes condition unless the PR description explicitly justifies the deviation. This protects QC-report completeness, run-metadata writes, and assertion placement against erosion.
- For multi-project work, PI uses git branches or project-specific directories — there is no framework support for "spawning a new analysis from a template", because that adds a workflow layer the framework deliberately does not own.

### Comment language: Chinese, with key "why" explanations

**All comments in both notebooks and src/ are in Chinese** (ADR-0009). The project's primary users are PI and non-CS students who read Chinese. Comments explain **why** a step is necessary, not just what it does:

- Why `target_sum=1e4` (library-size normalization standard)
- Why `seurat_v3` HVG (works on raw counts vs `seurat` on log-normalized)
- Why bidirectional gene ID sync is needed (scanpy wants symbols, cross-database queries need Ensembl)
- Why `compression="lzf"` (faster than gzip, ~30% smaller than uncompressed, preserves sparse CSR)
- Why scVI needs `>=20` epochs (variational inference needs enough iterations to converge)

Professional terms keep their English originals with Chinese explanations (e.g. "高可变基因 HVG", "批次效应 batch effect"). The goal: a student opening any notebook should understand every step within a few minutes of reading, without needing to consult external docs or an experienced bioinformatician.

---
