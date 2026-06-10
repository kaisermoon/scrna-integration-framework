# scRNA-seq Integration Framework

A modular scRNA-seq integration analysis framework for cross-disease, multi-source single-cell data analysis. Personal research infrastructure used by PI and students.

## What this framework is for

The framework's purpose is to enable PI and supervised students to produce publishable **biological findings** across disease projects (gastric precancerous lesions, RA, PCOS, VVC, ...). It is not a methods contribution — there is no software paper, no PyPI release, no plan to court external users. This framing is deliberate: it rules out the design pressure that comes with public-tool ambitions (broad API stability, exhaustive docs, generic extensibility).

The framework's code surface is intentionally tiny — two Python functions that fill genuine gaps in `scanpy` / `anndata` / `scverse`, plus directly-callable scorer functions in `scrna_integration.scorers`, plus conventions for everything else:

```python
from scrna_integration import read_with_manifest, load_markers
from scrna_integration.scorers import integration_metrics  # called directly in for loops
```

Anything else (run metadata, lineage, QC heterogeneity records, disease-ontology lookup, stage reports) is plain `adata.uns[...]` writes, plain `yaml.safe_load`, plain `sc.pl.*`, and notebook cells PI / students edit directly. There is **no `si.*` namespace tree**.

## Where to read more

- **`SPEC.md`** — implementation details: function signatures, manifest YAML schema, stage notebook cell sequences, R bridge idioms, memory discipline conventions, cross-method comparison protocol.
- **`docs/adr/`** — architectural decisions and their rationale (0001 thin framework over scanpy / 0002 rpy2 R bridge — *superseded by 0007* / 0003 plain code over plugin systems / 0004 framework deletion log + reviewer cheatsheet / 0005 load_markers as third function / 0006 YAML manifest for dataset facts / 0007 R bridge split by tool / 0008 absorb student-code by rewriting).
- **`_plan.md`** — current PR plan and project status.
- **`_memory.md`** — active project state, recent decisions, PI todos.
- **`README.md`** — installation, environment setup, running the GCPL pilot.

## Language

<!-- Domain glossary. One or two sentences max per term. Terms unique to this project only — general programming concepts don't belong. -->

**Stage**:
A coarse-grained pipeline phase (1–7, plus 6.5) that produces one h5ad checkpoint by file-naming convention. Stages are not framework objects; they are folders / file prefixes / stage notebooks.
_Avoid_: Step (reserved for finer intra-stage operations like "run PCA"), Phase (reserved for project lifecycle phases — planning/experiment/writing).

**Method run**:
A single invocation of a tunable method (e.g. one Harmony run, one scVI training). Multiple runs of the same step coexist in one h5ad as parallel `obsm` / `obs` slots, named by scanpy convention or PI's chosen suffix.
_Avoid_: Trial, attempt, experiment.

**Scorers**:
Directly-callable metric functions in `scrna_integration.scorers` (ADR-0009). For notebooks: `from scrna_integration.scorers import integration_metrics`. Called directly in explicit for loops — no callbacks. Replaces the removed `sweep()` function.
_Avoid_: Sweep function, grid search helper.

**Version** (of a stage h5ad):
A re-run of a stage that produces a new h5ad file with bumped version suffix (`_v2`, `_v3`). Old versions are never overwritten — file-system convention only.
_Avoid_: Revision, iteration.

**Source dataset**:
A single external data origin (one GEO accession, one published study, one cellranger output directory). Always recorded in `obs.source_dataset`. Multiple source datasets are integrated into one project.
_Avoid_: Dataset (ambiguous), Study.

**Project** (in framework sense):
A single integration analysis run end-to-end (e.g. "GCPL gastric precancerous integration"). One project consumes multiple source datasets and produces one set of stage h5ads. Distinct from Obsidian-vault project tagging.
_Avoid_: Analysis, pipeline run.

**Disease system**:
The broad research domain a project addresses (`"gastric"`, `"synovium"`, `"ovary"`, `"vaginal"`, ...). Required Layer 1 obs field; primary grouping dimension for cross-project analyses across PI's portfolio.
_Avoid_: Disease area, indication.

**Manifest**:
A per-source-dataset YAML file at `data/{source_dataset}/manifest.yaml` describing how to read the raw matrix, map its obs columns, join clinical metadata, and declare preprocessing state. Manifests are project assets; they live with the data and are version-controlled.
_Avoid_: Config, descriptor, recipe.

**Clinical metadata**:
External tables (xlsx/csv) supplying patient / sample-level fields not present in the matrix. Joined into obs by manifest declaration.
_Avoid_: Patient data, supplementary data.

**Original author annotation**:
Cell-type labels carried by a source dataset at ingest time. Stored in `cell_type_original_{source_dataset}_v1[_{role}]` columns, sparse across source datasets after integration.
_Avoid_: Pre-existing annotation, prior label.

**Cross-method comparison** (annotation; 中文：交叉比对):
The 06 notebook step where multiple annotation methods (default co-run: marker / mLLMCelltype / gene-set scoring; scANVI when a reference atlas exists; CellTypist as a commented-out candidate; plus original author annotations) are compared via confusion matrices, Sankey diagrams, and LLM verdicts to reach a final per-cluster cell-type assignment. Running several methods together is for cross-validation, not redundancy. Implemented in `06_annotated.ipynb`, not in framework code.
_Avoid_: Consensus (overloaded — LLM consensus is one input to the cross-method comparison), reconciliation, crosswalk (rejected as obscure jargon).

**Marker library**:
The `references/markers/` directory of CSV files. Long-term cross-project research asset. Loaded via `load_markers()` (see ADR-0005).
_Avoid_: Marker list, gene panel, signature.

**Disease ontology**:
A per-system YAML at `references/disease_ontology/{system}.yaml`. Loaded directly via `yaml.safe_load`.
_Avoid_: Disease tree, disease taxonomy.

**Promoted / Experimental / Deprecated**:
Convention values for `adata.uns["status"]`. PI writes the value directly when reviewing a re-run. There is no API for setting / querying these.
_Avoid_: Approved, blessed, retired.

**Layer 1 / Layer 2 / Layer 3** (obs schema):
The three tiers of obs columns. Layer 1 is required core (manifest enforces); Layer 2 is CellxGene-aligned with strong-warn + LLM best-effort fix; Layer 3 is project-defined free fields. See SPEC.md "obs Schema" section for the column lists.
_Avoid_: Required / recommended / optional (less precise — Layer 3 is "optional" only in the sense of not framework-prescribed; it can still be project-required).

**LLM verdict**:
The natural-language judgement an LLM writes for one cluster during 06 cross-method comparison, summarising what each annotation method called this cluster, what the markers / scores show, and what cell type the LLM recommends as final. PI reads every verdict and decides.
_Avoid_: AI annotation, LLM call (too generic).

**Sweep recommendations**:
The auto-generated section at the end of a 06 report listing parameter combinations the LLM judges worth running based on single-run results (e.g. clusters with high uncertainty, methods with systematic disagreement). PI decides whether to act. Note: named "sweep" as a workflow concept; the `sweep()` framework function was removed per ADR-0009.
_Avoid_: Suggestions, hints.

**Self-check**:
The one-line assertion on `adata.X` (sparse + float32) that every stage notebook places immediately before `adata.write_h5ad(...)` to catch dense / dtype regressions.
_Avoid_: Validation, sanity check (too generic).
