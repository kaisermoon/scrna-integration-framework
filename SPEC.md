# Implementation Specification

Implementation details for the scRNA-seq Integration Framework. Companion to `CONTEXT.md` (project glossary) and `docs/adr/` (architectural decisions). This document specifies *what is built and how it behaves*; CONTEXT.md defines *what terms mean*.

---

## Positioning

**The framework is a tool for biological discovery, not a methods contribution.** Its purpose is to enable PI and students to produce publishable biological findings across disease projects (gastric precancerous lesions, RA, PCOS, VVC, ...) — not to publish a software paper or release a community library. This framing rules out the design pressures that come with public-tool ambitions (broad API stability, exhaustive docs, generic extensibility).

**Stage 1 (current)**: Personal research infrastructure. Primary users are PI and supervised students. API may iterate aggressively; no public stability guarantees; lean docs. Function signatures may change between PRs during this stage — each signature change is logged in `_memory.md` so future readers can trace when an API shifted.

**Stage 2 (later)**: Re-evaluated when the framework has produced **at least one publishable biological finding from the GCPL pilot** (PR-4 complete + a stage 7 downstream module completes a real analysis). What "Stage 2" means at that point is decided then — it is not committed in advance to be a software paper, a public release, or any specific external deliverable.

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

## Architectural Stance: Three Functions + Conventions

**The framework code surface is intentionally tiny.** Three Python functions cover the genuine gaps in scanpy/anndata/scverse; everything else is convention plus reference data, used directly with `pandas` / `pyyaml` / `scanpy` / `anndata` native APIs.

```python
from scrna_integration import read_with_manifest, sweep, load_markers
```

These are the **only** three imports a notebook ever needs from the framework. Anything else (lineage, disease ontology lookup, QC skip records, run metadata, stage reports) is plain `adata.uns[...]` writes, plain `yaml.safe_load`, plain `sc.pl.*`, or notebook cells PI / students edit directly.

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

The framework's `__init__.py` re-exports `read_with_manifest`, `sweep`, and `load_markers`. There is no `si.*` sub-namespace tree.

### "Three functions" is the current reality, not an eternal promise

The current framework surface reflects today's understanding of where scanpy genuinely leaves a gap. A fourth (or further) framework function may be added when **all three** of the following are true:

1. **The same boilerplate appears in ≥3 stage notebooks** with non-trivial duplication (not just a 1-liner).
2. **The naive approach has been tried in real PR work and demonstrably caused maintenance burden** — copy-paste drift across notebooks, missed updates when the pattern needs to evolve, etc. "I anticipate this will become annoying" is not enough; the maintenance cost must already be observed (or — as with `load_markers` per ADR-0005 — PI confirms from prior real-world experience that the duplication is observed in lived practice, not anticipated).
3. **A new ADR is written** explaining why this case is an exception to ADR-0001 / 0003 / 0004 — what the naive approach was, why it failed, and what the new function does.

This bar is intentionally high. Most candidate "abstractions" that emerge during PR work fail criterion 2 — perceived future pain is not a license to add framework code. When in doubt, accept the duplication and revisit later.

If the framework grows from 3 functions to 7+, that is a signal the bar is being abused, not a signal the bar is wrong — tighten the bar (e.g. require ≥5 notebooks or PI sign-off in addition to the ADR).

---

## The Three Functions

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
11. **Records (does not load) raw matrix path** for SoupX: if manifest provides `input.raw_path`, writes it to `adata.uns["raw_matrix_path"]`. Stage 2 SoupX reads it on demand to avoid doubling stage 1 memory.
12. Computes baseline QC metrics (`n_genes`, `total_counts`, `pct_counts_mt`, `pct_counts_ribo`) on `adata.X` so they are present and column-aligned across all source datasets regardless of preprocessing state. These are stage-2 prerequisites; computing them at ingest avoids a missing-column branch at stage 2.
13. Returns plain `AnnData`. **Caller does anything they want next** — `adata.write_h5ad(...)` to checkpoint, scanpy ops to continue, etc.

The function is one file (`src/scrna_integration/io.py`) + minimal helpers. It does **not** call `validate_obs`, does **not** push metadata into a hidden namespace beyond `species` / `raw_matrix_path`, does **not** assign a "stage" tag. PI inspects obs after reading and decides whether the schema is OK.

### `sweep(fn, adata, candidates, scorer, output_dir) -> pd.DataFrame`

The genuine sweep gap: scanpy doesn't ship a generic "run this callable across a parameter grid, score each, write a comparison report" helper. This function:

1. Iterates over the Cartesian product of `candidates` (a dict of param-name → list of values).
2. Calls `fn(adata_copy, **params)` for each combo. (`fn` is any callable — `sc.pp.filter_cells`, a scvi-tools method, a user function — no wrapping required.)
3. Calls `scorer(adata_after, adata_before, params)` returning a dict of metrics.
4. Collates results into a `pd.DataFrame` (one row per combo, columns = params + metric keys).
5. Writes a markdown report to `output_dir` with the table and per-combo figures the scorer chose to save.
6. Returns the DataFrame. **PI sorts / filters / picks however they want.**

`scorer` is a plain function. The "composite scorer" idea (mandatory metrics + optional rank + optional LLM judge) becomes: PI either writes one scorer that returns whatever they need, or composes several scorers manually outside `sweep()`. The framework doesn't enforce shape.

### `load_markers(csv_path, roles=("canonical", "optional")) -> dict`

A marker-library loader for the `references/markers/*.csv` corpus. Justified per ADR-0005: PI's accumulated experience confirms the load + filter + groupby boilerplate appears in stage 6 annotation, per-cluster profiling, and downstream gene-set scoring notebooks; centralising the `role` semantics (`canonical` / `optional` / `negative`) prevents misuse.

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
stage 1    results/{project}_stage1_loaded_v{N}.h5ad     ingest + obs schema (read_with_manifest)
stage 2    results/{project}_stage2_qcd_v{N}.h5ad        QC + filter + doublet (scrublet) + ambient (SoupX)
stage 3    results/{project}_stage3_normalized_v{N}.h5ad normalize + log + HVG
stage 4    results/{project}_stage4_embedded_v{N}.h5ad   PCA / Harmony / scVI / scANVI (sweep + integration QC inline)
stage 5    results/{project}_stage5_clustered_v{N}.h5ad  multi-resolution Leiden
stage 6    results/{project}_stage6_annotated_v{N}.h5ad  multi-method annotation + cross-method comparison
stage 6.5  results/{project}_stage6_5_subset_v{N}.h5ad   subset re-cluster (re-run stage 3-5 on a subset; e.g. all T cells)
stage 7    results/downstream/{module}_v{N}.h5ad         downstream modules (one h5ad per module, see below)
```

### Stage 4 sweep includes integration QC

No embedding method is best on every dataset, so the workflow is to **run all candidates, then compare** — never to pick a default method up front. Candidates are peers: `X_pca` (baseline), `X_pca_harmony`, `X_scVI`, `X_scANVI` (when a labelled reference exists), cellxgene_census pretrained scVI (optional), and any future method PI adds (e.g. scPoli / scArches-family methods). Each writes its own `obsm` slot; none is privileged.

Two decision inputs, both first-class:

1. **Visual inspection (primary in PI's actual workflow)** — `sc.pl.umap` for every candidate embedding, coloured by `sample_id`, `batch`, and `cell_type` (original-author or marker-based). PI eyeballs how well each integrates batches without over-mixing biologically distinct cells. This visual read is the main basis for the decision, not a formality.
2. **Integration metrics** — batch entropy / iLISI / scIB suite (silhouette by batch / silhouette by celltype / kBET / graph connectivity / etc.), embedded in the stage 4 sweep report via `sweep(fn=embedding_method, ..., scorer=integration_metrics)`. The metric table corroborates or complicates the visual read.

PI weighs both and marks one embedding `promoted`. cellxgene_census pretrained scVI and scANVI are optional candidate cells (each with its prerequisite noted — census model coverage for the tissue / a labelled reference atlas respectively), not framework-prescribed steps.

**Adding an embedding method (extension pattern).** Same scanpy-native "parallel slots" mechanism as stage-5 clustering: a new method writes `adata.obsm["X_{method}"]` and is appended to the sweep's `use_rep` candidate list. No framework change — one added cell. This is how PI plugs in new integration methods as they are published ("学术追新").

### Stage 6.5 subset analysis

After stage 6 PI may want to refine a specific cell-type subset (all T cells, all epithelial cells, all SPEM-spectrum cells, ...). Stage 6.5 re-runs stage 3 (HVG re-selection on the subset) → stage 4 (re-embedding) → stage 5 (re-clustering) → stage 6 (re-annotation) on the subset. Output goes to `stage6_5_subset_v{N}.h5ad` with `adata.uns["subset_of"] = "results/gcpl_stage6_annotated_v1.h5ad"` and `adata.uns["subset_filter"] = "cell_type_final_v1.isin(['CD4 T', 'CD8 T', 'Treg'])"` for traceability.

Stage 6.5 is implemented as a single notebook (`stage6_5_subset.ipynb`) that calls the same scanpy / framework API as stages 3-6 on the subset. Re-clustering parameters are project-specific (the T-cell subset usually needs different HVG / resolution than the global UMAP).

### Stage 7 downstream modules

Stage 7 is a fan-out — multiple downstream modules consume stage 6 (or stage 6.5) output independently and run in parallel. Each module has its own h5ad and its own notebook:

| Module | Notebook | Tool stack | Status | Student-code reference (re-implement per ADR-0008) |
|---|---|---|---|---|
| Differential expression (per-cluster, scanpy native) | `notebooks/stage7/deg.ipynb` | `sc.tl.rank_genes_groups` | **PR-3 in scope** | legacy-GCPL `08_differential_expression.ipynb` |
| Pseudobulk DEG (cross-condition: disease vs control) | `notebooks/stage7/pseudobulk_deg.ipynb` | DESeq2 (subprocess Rscript) | **PR-3 in scope** | — (decoupler pseudobulk + DESeq2 contrast) |
| CNV inference | `notebooks/stage7/cnv.ipynb` | infercnvpy (pure Python) | **PR-3 in scope** | `workflow_for_pseudotime/4.2_*` (gene positions) + `4.3_*` (infercnvpy run) |
| Pathway enrichment | `notebooks/stage7/pathway.ipynb` | GSEApy / decoupler / Reactome | PR-5+ | CD4 deep-analysis template GSEA part (reference only) |
| Pseudotime + root identification | `notebooks/stage7/pseudotime.ipynb` | CytoTRACE (cellrank) / Monocle3 (Rscript) / transcriptome entropy (numpy) | PR-5+ | `4.3_*` (entropy/CytoTRACE/Monocle3 export) + `4.4_*`/`4.5_*` (multi-metric root id) + `11.2_*` (lineage Monocle3) |
| GRN | `notebooks/stage7/grn.ipynb` | pySCENIC | PR-5+ | — |
| Cell communication | `notebooks/stage7/cell_communication.ipynb` | CellChat / NicheNet / CellPhoneDB | PR-5+ | CD4 deep-analysis LR part (reference only) |
| Differential abundance | `notebooks/stage7/abundance.ipynb` | scCODA / Milo | PR-5+ | `11_all_celltype_proportion_analyse.ipynb` (scCODA + Mann-Whitney + Cliff's delta + effect sizes) |
| Gene co-expression modules | `notebooks/stage7/gene_modules.ipynb` | hdWGCNA (subprocess Rscript) | PR-5+ | — |

The 3 modules in PR-3 scope (DEG / pseudobulk DEG / CNV) form the **minimum set needed to reach the GCPL pilot's first biological finding** — distinguishing tumor vs normal cells, finding genes differential along the CAG → IM → dysplasia axis, and identifying disease-vs-control pseudobulk DEGs. The remaining 6 modules are added in subsequent PRs as the GCPL analysis surfaces specific scientific questions.

Multiple method runs at any stage coexist in the same h5ad as additional `obsm` / `obs` / `var` slots, following scanpy's natural anndata usage. Re-runs bump the version suffix; old versions remain on disk for rollback.

### Lineage = file naming + uns annotation

A re-run writes its own provenance directly into `adata.uns`:

```python
adata.uns["status"] = "promoted"          # or "experimental" / "deprecated"
adata.uns["upstream"] = ["results/gcpl_stage3_normalized_v1.h5ad"]
adata.uns["notes"] = "tried scVI in place of Harmony; not promoted"
```

PI inspects these dict entries directly. There is no API for `promote()` / `deprecate()` / `show_dependents()`. When PI wants to know "what depends on stage4_v1", they `glob` over `results/` and read each file's `adata.uns["upstream"]`.

### Default upstream selection

There is no automatic default. Every notebook has a parameters cell where the upstream path is **explicitly named**:

```python
# === PARAMS ===
upstream_path = "results/gcpl_stage3_normalized_v1.h5ad"
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

**Example 2 — recording a Harmony embedding run** at stage 4:

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

**Example 3 — recording PI's `cell_type_final` decision** at stage 6:

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
- `tissue` / `tissue_ontology_term_id` — tissue-aware marker selection, stage 6 annotation, cross-tissue comparison
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
                                      # for stage 2 SoupX. Stage 1 does NOT load it (avoids 2x memory).
                                      # Absent → SoupX skipped at stage 2.
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

## R Bridge: split by tool

R interoperability is unavoidable, but the mechanism is chosen **per tool** — see ADR-0007. There is no single canonical path and no `_r_bridge` framework module. Three cases:

1. **Pure Python, no R** — InferCNV (`infercnvpy`: `cnv.tl.infercnv` / `cnv.tl.cnv_score`) and CytoTRACE (`cellrank` `CytoTRACEKernel`). These have mature Python packages validated in `student-code`; do not reach for R.
2. **rpy2 `%%R` in-notebook** — SoupX at stage 2. Per-sample matrices are small and conversion is stable; keeps QC coherent inside one notebook. Same style as `legacy-GCPL/02_quality_control.ipynb`.
3. **subprocess `Rscript` + temp files** — Monocle3, UCell, DESeq2, hdWGCNA, and other heavy / R-only tools. Write `.mtx`/`.csv`, call `subprocess.run(["Rscript", "--vanilla", ...])`, read results back. This is the **primary** path for heavy R tools because rpy2 + anndata2ri is brittle on large-object conversion and breaks on dependency upgrades; the students' working downstream code uses subprocess throughout.

The only framework-side R concession is that `read_with_manifest` handles `format: "rds"` internally (otherwise every IO call site needs rpy2 ceremony just to read a Seurat object). All other R access follows the split above, written directly in notebooks.

### rpy2 idiom — SoupX (stage 2, small in-notebook exchange)

Three lines at the top of the QC notebook:

```python
%load_ext rpy2.ipython
import anndata2ri
anndata2ri.set_ipython_converter()
from rpy2.robjects.packages import importr
```

After these, `%R` and `%%R` cells work and AnnData ↔ SingleCellExperiment conversion happens transparently in `localconverter` blocks.

### Concrete example — SoupX (stage 2 ambient correction)

```python
# Cell A: load the SoupX library and read raw matrix on demand
importr('SoupX')
adata_raw = sc.read_10x_mtx(adata.uns["raw_matrix_path"])   # not loaded at stage 1; loaded here

# Cell B: run SoupX per sample (R cell — define and call the R function)
%%R -i counts_matrix -i raw_counts_matrix -i cluster_labels -o corrected_counts_matrix
library(SoupX)
sc <- SoupChannel(raw_counts_matrix, counts_matrix)
sc <- setClusters(sc, cluster_labels)
sc <- autoEstCont(sc)
corrected_counts_matrix <- adjustCounts(sc)
```

PI loops over samples in Python, passes each sample's matrices in via `-i`, gets corrected counts back via `-o`. Same pattern as `legacy-GCPL/02_quality_control.ipynb`.

### Concrete example — Monocle3 (stage 7 pseudotime, subprocess Rscript)

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

### Concrete example — pseudobulk DESeq2 (stage 7 cross-condition DEG, subprocess Rscript)

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

**rpy2 (SoupX only):**
- **Always wrap conversion in `localconverter`** — anndata2ri's converter is not always active by default; the explicit context manager prevents surprise type errors.
- **Use `%%R -i input1 -i input2 -o output`** to declare exactly which Python variables enter and leave the R cell. No implicit globals across the language boundary.
- **Never `library(...)` inside a hot loop** — `importr(...)` once at the top.

**subprocess (heavy R tools):**
- **The `.R` script lives in `scripts/`**, takes a working-directory path (and any contrast args) on the command line, and is independently runnable — debug it with `Rscript --vanilla scripts/foo.R results/_tmp` outside Jupyter.
- **Always pass `check=True`** to `subprocess.run` so a failing R script raises in Python instead of silently continuing.
- **Read R-produced figures back into the notebook** (`IPython.display.Image`) so results stay visible in the notebook even though computation happened out of process.
- **Temp dirs go under `results/_<tool>_tmp/`** (gitignored) and are cleaned at the end of the cell when large.

### When rpy2 conversion breaks

If SoupX's rpy2 path fails after a numpy/scipy/anndata upgrade, fall back to the subprocess pattern for SoupX too (export per-sample matrices, call an `.R` script). This is the per-incident workaround for the one rpy2 site; the heavy-R tools are already on subprocess by default.

### Environment management

Two conda environment files at repo root: `environment.yml` (Python) and `environment-r.yml` (R). The R env is structured by stage so users skip optional packages (e.g. someone only doing QC skips Monocle3). README ships a one-liner `mamba env create + update` for both.

**环境隔离硬规定（agent 与 PI 都必须遵守）**：所有包安装一律在**专用命名 conda 环境**内进行，**绝不允许动 base / 主环境**。

- Python 环境固定命名 `scrna-integration`，R 环境固定命名 `scrna-integration-r`（写入 `environment.yml` / `environment-r.yml` 的 `name:` 字段）。
- agent **禁止** `pip install` / `conda install` 到 base 或任意非专用环境；禁止 `conda install -n base`；禁止无隔离的全局 `pip install`。
- 任何安装前必须先 `conda activate scrna-integration`（或 `-r`）；环境不存在则先 `conda env create -f environment.yml`。
- CI 不实际装重型科学栈——`pyproject.toml` 与两个 `environment.yml` 是**声明**，真实安装由 PI/coder 在本地专用环境完成。验收以"依赖声明完整、smoke test（纯 import 框架包）通过、环境文件可被 conda 解析"为准，不要求 CI 跑完整安装。
- reviewer 审 env 相关 PR 时红线检查：diff 出现写向 base/全局/非专用环境的安装命令或文档指引 → flag。

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
| `references/cell_type_ontology/{system}.yaml` (cell type hierarchy) | LLM in stage 6 cross-method comparison already handles `"T cell"` vs `"CD8+ Tem"` semantic alignment without an explicit ontology. Maintaining a parallel cell type hierarchy adds cost without observed benefit. |
| `references/genesets/` (gene set library) | GSEApy / decoupler ship MSigDB / KEGG / Reactome built-in. Project-specific gene sets (e.g. SPEM transition signature) are best stored where they're used — in the relevant analysis notebook or per-project data dir — not promoted to a framework-wide library. |
| `references/embeddings/` (cellxgene_census ref metadata) | Stage 4 notebooks call `cellxgene_census` API directly with hard-coded metadata. Locally archiving model versions adds maintenance without solving a real problem. |
| `references/prompts/` (LLM verdict / annotation prompts) | Per-stage LLM prompts live in the stage notebook scripts directly. Centralising them creates a sync burden between notebook and prompt file with no clear advantage. |

The decision rule for adding a third asset category to `references/` mirrors ADR-0001 / 0003 / 0004: a candidate must demonstrate observed real-world benefit (not anticipated convenience) and a new ADR before being added.

---

## QC Heterogeneity

Some source datasets ship pre-filtered (Tsubosaka, Nowicki-Osuch, Kim) — re-running scrublet / SoupX on already-cleaned data is a scientific error. Manifest declares preprocessing state via `preprocessing_done` and per-step `qc_overrides`. PI's QC notebook then **skips** matching steps based on the manifest dict — the notebook reads the manifest, branches on each step, runs scanpy QC functions where applicable.

```python
# In stage 2 QC notebook
manifest = yaml.safe_load(open(f"data/{source_dataset}/manifest.yaml"))
preprocessing_done = manifest.get("preprocessing_done", [])
qc_overrides = manifest.get("qc_overrides", {})

if "doublet_removal" not in preprocessing_done and not qc_overrides.get("doublet_removal", {}).get("skip"):
    sc.external.pp.scrublet(adata)
else:
    print(f"Skipping doublet_removal: {qc_overrides.get('doublet_removal', {}).get('reason', 'already done by author')}")
    adata.obs["doublet_score"] = float("nan")    # column-align across datasets
    adata.obs["predicted_doublet"] = False
```

The framework does **not** ship a `qc_runner` function. The stage 2 notebook (`stage2_qcd.ipynb`) contains the branch logic above directly; PI edits PARAMS at the top and runs.

### Column alignment

For final cross-dataset analyses to work, QC-derived obs columns must be present and aligned across all integrated cells:

- `obs.n_genes`, `obs.total_counts`, `obs.pct_counts_mt`, `obs.pct_counts_ribo` — always recomputed from `adata.X` regardless of preprocessing state. The stage 2 notebook includes these computations.
- `obs.doublet_score`, `obs.predicted_doublet` — present as columns; **NaN** on cells where doublet detection was skipped (encodes "not applicable", not "missing data").
- `obs.ambient_correction_applied` (bool) — true on cells where ambient correction ran; false elsewhere.

PI reads these conventions in the stage notebook and follows them.

### What the framework does NOT enforce

PI explicitly accepts that QC heterogeneity propagates into stage 4-7 results — strict cross-dataset alignment was considered and rejected:

- **Doublet alignment is NOT mandatory across datasets**. A dataset that ships already-doublet-removed (Tsubosaka, Nowicki, Kim) skips scrublet at stage 2; a raw cellranger dataset (Nancang) runs scrublet. The two datasets enter stage 3+ with different doublet histories. The framework does not force a "re-run scrublet on already-cleaned data to produce a comparable score" alignment step. PI accepts the resulting risk that Harmony/scVI batch correction at stage 4 may treat residual doublets in some datasets as `false batch effect`. Mitigation is a disclaimer in the stage 4 sweep report when `qc_heterogeneous=True`, not an alignment step.
- **Ambient correction is physically gated by data availability.** SoupX requires `input.raw_path` in the manifest (the cellranger raw_feature_bc_matrix). Datasets that ship only filtered matrices cannot run SoupX — there is no way to "force align". `obs.ambient_correction_applied` records the actual situation per dataset; downstream consumers read it.
- **Stage 7 cross-condition DEG / pseudobulk** carries the same disclaimer. Pseudobulk averaging dilutes single-cell-level QC heterogeneity to some extent (multi-cell mean per sample), but cannot eliminate systematic bias when one cohort has been doublet-removed and another has not. PI accepts this risk for now; if it produces visibly biased DEG results in the GCPL pilot, alignment can be reconsidered as a per-project step.

### Cross-method comparison reads QC context

When stage 6 cross-method comparison includes original author annotations (`cell_type_original_{source_dataset}_v1`), the LLM verdict prompt is supplied with each source dataset's `qc_skipped` record so the LLM can interpret author labels in their proper QC context. For example, when comparing `cell_type_original_Tsubosaka_2023_v1 = "T cell"` (annotated on a doublet-removed matrix) against `cell_type_llm_v1 = "T cell + contaminated"` (annotated on a matrix that still contains doublets), the LLM should recognise the difference may be a QC artefact rather than a true biological disagreement. This is implemented in the `stage6_annotated.ipynb` LLM verdict prompt construction — not in framework code.

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

## Stage 6: Multi-method Annotation in the Notebook

Stage 6 is the most distinctive scientific layer of the framework, but it is implemented as a **notebook** (`stage6_annotated.ipynb`), not a framework module. Annotation methods run as **independent cells in the notebook** — PI is free to choose any execution order; the framework does not anchor "correct order" between manual marker annotation and the automatic methods. (Order experimentation is encouraged: running automatic methods first lets PI cross-check intuition, while running manual marker first preserves PI's prior judgement free of automatic-method anchoring.)

**Running several methods in parallel is the point, not redundancy.** The most accurate per-cluster cell type comes from cross-comparing independent methods — agreement raises confidence, disagreement flags clusters needing PI attention. So the default is to run multiple methods together, not to pick one.

Methods fall into three tiers:

**Default co-run (run together every time):**

1. **Marker dotplot** — PI views dotplots of canonical/optional markers from `references/markers/{tissue}_{purpose}.csv` and manually labels each leiden cluster, writing `obs.cell_type_marker_v1`.
2. **mLLMCelltype consensus via OpenRouter** — writes `obs.cell_type_llm_v1` + `_uncertainty` + `_consensus_log`. Same code pattern as `legacy-GCPL/06_annotation.ipynb`.
3. **Gene-set scoring** — AUCell / UCell / scanpy `score_genes` over canonical marker sets, producing continuous scores in `obs.score_*` columns. Cross-cluster summarisation (mean score + pct-positive per cluster, transition/mixed flagging) follows the pattern in `student-code/.../6.3_UCell_cluster_mean_score&pct_pos.py`. Used as supplementary evidence in cross-method comparison rather than as standalone cell-type labels.

**Default-active when a reference atlas is available (data-availability gated, not preference):**

4. **scANVI label transfer** — writes `obs.cell_type_scanvi_v1` + `_uncertainty` from a labelled reference atlas. The cell is active by default; it runs whenever a suitable reference atlas exists for the disease system. When none is available the cell is skipped (it is the only default method with a hard external prerequisite — a labelled atlas for this tissue). The notebook cell header states the prerequisite explicitly.

**Candidate (commented-out cell, enable when applicable):**

5. **CellTypist pretrained classifier** — writes `obs.cell_type_celltypist_v1` + `_uncertainty`. Left as a commented-out cell with a header noting "enable when a CellTypist model exists for this tissue." Not run by default because PI's gastric / synovial systems do not currently have well-matched pretrained models; enabling it is a one-line uncomment when a relevant model appears.

After the method cells run, the notebook constructs the cross-method comparison over whichever methods produced labels:

6. **Confusion matrices + Sankey diagrams** between method pairs (with Cohen's kappa) using scanpy / seaborn natively.
7. **LLM verdict per cluster** — the notebook constructs a prompt containing all available method labels (including original-author annotations from `cell_type_original_*` if present), top markers from `sc.tl.rank_genes_groups`, gene-set score profile, marker_db coverage, and `adata.uns["qc_skipped"]` per source dataset. OpenRouter is called directly. The verdict for each cluster is written to a per-cluster markdown file.
8. **PI reads every cluster's verdict and decides** — `obs.cell_type_final_v1` is set manually. With one integrated adata per project (per-project N is typically 15–30 leiden clusters), reading every verdict takes 30–60 minutes and is the standard workflow. There is no auto-pick by confidence threshold — every cluster's final label is PI's call.

**There is no `si.report.stage6_annotation` function.** The notebook is the unit of work. When new annotation methods appear, PI / students add a cell to `stage6_annotated.ipynb` and rename the obs column following `cell_type_{method}_v{N}` naming.

### `cell_type_final_v{N}` versioning

PI may revisit annotations weeks later (e.g. after seeing a stage 7 result that suggests refining cluster 7 from `SPEM` into `SPEM_complete` / `SPEM_incomplete`). Following the same convention as stage h5ad versioning, **the previous `cell_type_final_v{N}` column is never overwritten** — a new `cell_type_final_v2` column is added with the revised labels, and `adata.uns["cell_type_final_v2_notes"]` records the revision rationale. PI / agents reference whichever version is appropriate for the analysis at hand.

### Stage 6.5 subset re-annotation reflows back

When stage 6.5 produces refined sub-cluster annotations (e.g. CD4/CD8/Treg/MAIT from a T-cell subset re-clustering), the refined labels are reflowed into the main stage 6 adata as a new column `cell_type_final_subset_v1`. Cells in the refined subset receive the fine label; cells outside the subset (epithelial, stromal, etc.) receive `NaN` in this column. The original `cell_type_final_v1` column remains untouched. This preserves the granularity hierarchy (broad in `_final_v1`, fine in `_final_subset_v1`) without polluting either, and downstream stages can choose the column appropriate to their question.

### Original author annotations

Manifest's `original_annotations` section causes IO to rename author labels to `cell_type_original_{source_dataset}_v1[_{role}]`. After multi-source integration these columns are sparse (NaN on cells from datasets without that annotation). The stage 6 notebook treats them as additional methods in the cross-method comparison — the LLM verdict prompt naturally references them when present and reads `adata.uns["qc_skipped"]` so the LLM can interpret author labels in their proper QC context (see "Cross-method comparison reads QC context" in the QC Heterogeneity section).

### Sweep recommendations

The stage 6 notebook ends with one final cell that constructs an LLM call to generate **sweep recommendations** based on the single-run results: which clusters had high uncertainty across methods, which methods showed systematic disagreement, which marker_db gaps the LLM observed. Output is appended to the stage 6 markdown report. PI reads recommendations and decides whether to invoke `sweep()` on annotation parameters (e.g. different LLM model sets, different scANVI references). Cost is one OpenRouter call per stage 6 run — acceptable. This is the operational form of the project's "teach AI to make initial scRNA-seq judgements" thesis.

### Per-cluster deep profile

A second notebook (`stage6_per_cluster.ipynb`) loops over clusters and emits one markdown per cluster with UMAP highlight, top-marker dotplot, gene-set score violins, cross-disease abundance boxplot, and an LLM-written narrative paragraph. **No framework class, no panel registry — just a notebook with a `for cluster in adata.obs[label_col].cat.categories: ...` loop.**

---

## Notebooks (Directly Runnable, Not Templates)

The `notebooks/` directory contains **directly runnable notebooks**, not templates. PI edits the PARAMS cell at the top and runs all cells — no copy / rename ceremony. The framework's standardisation ships as the actual notebook content, not as a template-instantiation system.

```
notebooks/
├── stage1_loaded.ipynb              # PARAMS cell + read_with_manifest + write h5ad
├── stage2_qcd.ipynb                 # QC with manifest-driven skip logic + scrublet + SoupX (rpy2)
├── stage3_normalized.ipynb          # normalize + log + HVG
├── stage4_embedded.ipynb            # PCA + Harmony + scVI + scANVI + sweep with integration_metrics
├── stage5_clustered.ipynb           # multi-resolution Leiden + sweep
├── stage6_annotated.ipynb           # 5-method annotation + cross-method comparison + LLM verdict
├── stage6_per_cluster.ipynb         # per-cluster deep profile
├── stage6_5_subset.ipynb            # subset re-cluster (T cells / epithelial / etc.)
└── stage7/                          # downstream modules; each consumes stage 6 (or 6.5) output
    ├── deg.ipynb                    # PR-3 — sc.tl.rank_genes_groups
    ├── pseudobulk_deg.ipynb         # PR-3 — DESeq2 (rpy2)
    ├── cnv.ipynb                    # PR-3 — InferCNV (rpy2)
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
2. **Imports cell** (scanpy / framework imports / rpy2 if needed).
3. **Load upstream h5ad** (or call `read_with_manifest` for stage 1).
4. **Stage logic cells** — scanpy native APIs.
5. **Run-metadata cells** — `adata.uns[...] = {...}` writes documenting what was done.
6. **Report cells** — `sc.pl.*` + `plt.savefig` + per-stage markdown summary.
7. **Final cell** — `adata.write_h5ad(out_path, compression="lzf")` + `del adata; gc.collect()`.

### Per-notebook cell sequence (specification for PR-3 implementation)

The cell sequences below are the **specification** PR-3 coder agents implement against. Each line is one notebook cell; markdown cells are prefixed `[md]`, code cells are `[code]`.

#### `stage1_loaded.ipynb`

```
[md]   # Stage 1: Multi-source ingest + obs schema standardisation
[code] # === PARAMS ===
       MANIFEST_PATH = "data/Nancang_2025/manifest.yaml"
       OUTPUT_PATH   = "results/{project}_stage1_loaded_v1.h5ad"
       RANDOM_SEED   = 42
[code] # imports
[code] # call read_with_manifest, inspect returned adata
[code] # print obs.head() / var.head() / .uns keys for quick visual check
[code] # baseline QC metrics already on obs (read_with_manifest computed them); print summary
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write h5ad with compression="lzf", del adata, gc.collect
```

#### `stage2_qcd.ipynb`

```
[md]   # Stage 2: QC + filter + doublet (scrublet) + ambient (SoupX, conditional)
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, manifest path lookup
[code] # imports (scanpy + rpy2 if any source uses SoupX)
[code] # load adata + load manifest dict (for preprocessing_done / qc_overrides)
[md]   ## QC visualisation pre-filter
[code] # violins / scatters of n_genes / pct_mt / total_counts per source_dataset
[md]   ## Doublet detection (skipped per-dataset based on manifest)
[code] # for each source_dataset: if "doublet_removal" in preprocessing_done OR qc_overrides[doublet_removal].skip:
       #     fill obs.doublet_score = NaN, predicted_doublet = False, log skip reason
       # else:
       #     sc.external.pp.scrublet(adata_per_source); fold back into main adata
[md]   ## Ambient correction via SoupX (conditional on input.raw_path)
[code] # for each source_dataset: if adata.uns["raw_matrix_path"][source] is not None:
       #     load raw matrix, call SoupX (R cell), write corrected counts back
       #     obs.ambient_correction_applied = True
       # else:
       #     obs.ambient_correction_applied = False; log "no raw matrix"
[md]   ## Filter
[code] # apply min_genes / max_genes / max_pct_mt thresholds
[md]   ## QC visualisation post-filter (for before/after comparison)
[code] # same plots as pre-filter; PI eyeballs the change
[code] # adata.uns["qc_skipped"] structured record (per source_dataset, per step)
       # adata.uns["qc_heterogeneous"] = True if any skipped
       # adata.uns["filter_v1"] = {params, cells_in, cells_out, timestamp}
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write h5ad, del adata, gc.collect
```

#### `stage3_normalized.ipynb`

```
[md]   # Stage 3: normalize + log + HVG
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

#### `stage4_embedded.ipynb`

```
[md]   # Stage 4: multi-method embedding (all peers) + integration sweep
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
[md]   ## Sweep across embeddings with integration_metrics (corroborating metrics)
[code] # sweep(fn=neighbors_then_umap, adata=adata,
       #       candidates={"use_rep": [m for m in embedding_methods present]},
       #       scorer=integration_metrics, output_dir=...)
[code] # PI weighs UMAP visuals + metric table, decides which to mark promoted
[code] # adata.uns["harmony_v1"] / ["scvi_v1"] / etc. metadata writes
       # adata.uns["status"] = "experimental" (PI manually changes to "promoted" later)
[code] # cast all obsm latent matrices to float32
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write, del, gc
```

**Adding an embedding method**: write `adata.obsm["X_{method}"]` and add `"X_{method}"` to the sweep `use_rep` candidates — one cell, no framework change (see "Stage 4 sweep includes integration QC" above).
[code] # adata.uns["harmony_v1"] / ["scvi_v1"] / ["scanvi_v1"] metadata writes
       # adata.uns["status"] = "experimental" (PI manually changes to "promoted" later)
[code] # cast all obsm latent matrices to float32
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write, del, gc
```

#### `stage5_clustered.ipynb`

```
[md]   # Stage 5: multi-resolution Leiden + sweep
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, resolutions (list), use_rep
[code] # imports + load
[md]   ## Multi-resolution Leiden
[code] # for res in resolutions: sc.tl.leiden(adata, resolution=res, key_added=f"leiden_res_{res}")
[md]   ## Sweep across resolutions
[code] # sweep(fn=lambda adata, resolution: sc.tl.leiden(adata, resolution=resolution),
       #       candidates={"resolution": resolutions},
       #       scorer=clustering_metrics)  # silhouette + ARI per resolution
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

No framework code is needed for this. The extension mechanism is the same scanpy-native "multiple results coexist as parallel slots" pattern as stage-4 embeddings:

- Any clustering method writes its labels to its own obs column: `obs["acdc_clusters"]`, `obs["{method}_clusters"]`, etc. Multi-resolution Leiden already does this (`leiden_res_0.5`, `leiden_res_1.0`, ...).
- Methods that sweep a parameter (Leiden over resolution) go through `sweep()`. Methods that auto-pick (ACDC searches for an optimal partition itself) just write their single result column directly — no sweep needed.
- All resulting cluster columns are compared the same way: `sc.pl.umap` coloured by each, plus ARI / silhouette across columns. PI picks which column feeds stage 6.
- A new method is one added cell, not a framework change. ACDC specifically was tried in GCPL but never completed a run (prohibitively slow on PI's data), so it is **not** a default and `acdc_py` is **not** a PR-0a dependency; it is one example of the kind of method this slot accepts when PI chooses to install it.

#### `stage6_annotated.ipynb`

```
[md]   # Stage 6: Multi-method annotation + cross-method comparison
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, MARKER_CSV, leiden_col, OPENROUTER_API_KEY
[code] # imports + load + load_markers(MARKER_CSV)
[md]   ## Method 1: marker dotplot (PI manually labels)
[code] # sc.pl.dotplot(adata, var_names=markers, groupby=leiden_col)
[code] # PI inspects, fills marker_assignments dict, then:
       # adata.obs["cell_type_marker_v1"] = adata.obs[leiden_col].map(marker_assignments)
[md]   ## Method 2: mLLMCelltype consensus via OpenRouter
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
       # call OpenRouter, write per-cluster verdicts to results/figures/stage6_verdicts/
[md]   ## PI sets cell_type_final_v1 (manual review against verdicts)
[code] # pi_decisions = {leiden_id: cell_type, ...}
       # adata.obs["cell_type_final_v1"] = adata.obs[leiden_col].map(pi_decisions)
       # adata.uns["cell_type_final_v1_notes"] = {...}
[md]   ## Sweep recommendations (LLM-generated)
[code] # one OpenRouter call producing recommendations markdown appended to stage 6 report
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write, del, gc
```

#### `stage6_per_cluster.ipynb`

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

#### `stage6_5_subset.ipynb`

```
[md]   # Stage 6.5: subset re-cluster (e.g. T cells, SPEM-spectrum, epithelial)
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, subset_filter (str expression on obs), n_top_genes, resolutions
[code] # imports + load
[code] # adata_sub = adata[adata.obs.eval(subset_filter)].copy()
       # adata_sub.uns["subset_of"] = UPSTREAM_PATH
       # adata_sub.uns["subset_filter"] = subset_filter
[code] # re-run stage 3 (HVG re-selection on subset)
[code] # re-run stage 4 (re-embedding)
[code] # re-run stage 5 (re-clustering)
[code] # re-run stage 6 (re-annotation; usually finer cell types)
[code] # reflow refined labels back to main adata as cell_type_final_subset_v1 column
       # main_adata = sc.read_h5ad(UPSTREAM_PATH)
       # main_adata.obs["cell_type_final_subset_v1"] = ...  # fill subset cells, NaN elsewhere
       # main_adata.write_h5ad(UPSTREAM_PATH_with_subset_v1)  # bumped version of main h5ad
[code] # write subset h5ad + write updated main h5ad
```

#### `stage7/deg.ipynb` (PR-3 in scope)

```
[md]   # Stage 7: per-cluster differential expression
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, label_col, method
[code] # imports + load
[code] # sc.tl.rank_genes_groups(adata, groupby=label_col, method="wilcoxon")
[code] # sc.pl.rank_genes_groups + sc.pl.rank_genes_groups_dotplot
[code] # export top-N DEG per cluster as CSV
[code] # adata.uns["rank_genes_groups_v1"] is set by scanpy; record additional metadata
[code] # memory self-check: assert sp.issparse(adata.X) and adata.X.dtype == np.float32
[code] # write, del, gc
```

#### `stage7/pseudobulk_deg.ipynb` (PR-3 in scope)

```
[md]   # Stage 7: pseudobulk cross-condition DEG via DESeq2 (rpy2)
[code] # === PARAMS === UPSTREAM_PATH, OUTPUT_PATH, sample_col, contrast_col, contrast_levels
[code] # imports (decoupler for pseudobulk, rpy2 for DESeq2)
[code] # pseudobulk per (sample_id × cell_type_final_v1) using decoupler
[code] # for each cell_type: hand pseudobulk count matrix + sample metadata to DESeq2 (R cell)
[code] # collate results across cell types into one CSV / dict in adata.uns
[code] # plot: heatmap of top DEGs per cell type, volcano per contrast
[code] # write CSVs to results/, optionally adata.uns["pseudobulk_deg_v1"]
```

#### `stage7/cnv.ipynb` (PR-3 in scope)

```
[md]   # Stage 7: CNV inference via infercnvpy (pure Python, no R)
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
- Every notebook uses scanpy native APIs and the three framework functions (`read_with_manifest` / `sweep` / `load_markers`) only where they fill a real gap.
- Run-metadata writes follow the "Run Metadata" section above (plain `adata.uns[key] = {...}` with versioned keys).
- **Cell-sequence specs are binding.** When a notebook drifts from its spec (cells removed / reordered / renamed) during PR review, that is a request-changes condition unless the PR description explicitly justifies the deviation. This protects QC-report completeness, run-metadata writes, and assertion placement against erosion.
- For multi-project work, PI uses git branches or project-specific directories — there is no framework support for "spawning a new analysis from a template", because that adds a workflow layer the framework deliberately does not own.

---
