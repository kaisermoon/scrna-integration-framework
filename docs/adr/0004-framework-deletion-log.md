# Framework deletion log: from 11 si.* APIs to 2 functions

A second-round PI review during the grilling session triggered a major cut to the framework's code surface. Earlier iterations of `CONTEXT.md` introduced a `si.tracking` / `si.lineage` / `si.markers` / `si.disease` / `si.report` / `si.scorers` / `si.qc` / `si.memory` namespace tree — eleven public APIs in total. PI rejected every layer that wrapped existing scanpy / anndata / pandas / yaml functionality, on the grounds that "students learn scanpy first; any framework-side wrapping creates a parallel vocabulary they must learn on top". Inspection of `references/legacy-GCPL/` confirmed that PI's own real-world notebooks are entirely scanpy-native — there was no precedent in PI's codebase for the wrapper APIs the framework was proposing.

The framework was rewritten down to three functions (`read_with_manifest`, `sweep`, `load_markers`) plus conventions (file naming, manifest schema, `adata.uns` writes by PI directly, directly runnable per-stage notebooks, reference data CSV/YAML loaded via `load_markers` and stock `pyyaml`). The third function `load_markers` was added later under ADR-0005 as a documented exception to the bar.

## Considered Options

- **Keep the wrapper APIs but document them better**: rejected. The cost was readability and learning surface, not documentation. Better docs do not address the core complaint.
- **Keep wrappers for "common patterns" only** (`si.markers.load`, `si.disease.ancestors`): rejected. Each wrapper is one more API to remember, one more breaking-change surface during framework iteration, and one more place where PI / agents wonder "do I use the wrapper or the underlying function". The savings (one line of pandas / yaml code) are not worth the cognitive load.
- **Keep `si.lineage` because lineage tracking is novel**: rejected. The functionality reduces to (a) writing `adata.uns["status"]`, (b) writing `adata.uns["upstream"]`, (c) using `glob` over `results/`. None of the three benefits from being wrapped.

## Consequences

- The framework's code surface is two Python files (`io.py` + `sweep.py`) plus reference data. New contributors / agents read the entire framework in under an hour.
- Conventions that were going to be enforced by runtime checks (`validate_obs`, `memory.check`, `lineage.show_dependents`) are now enforced by code review (the `code-reviewer` agent reads `CONTEXT.md` and flags violations on PR diffs).
- Stage notebooks (`notebooks/stage1_loaded.ipynb` … `stage6_per_cluster.ipynb`) carry significantly more weight in the framework's design — they're now the primary unit through which PI standardises analyses across projects. PI's later third-round feedback explicitly rejected the `_template_` prefix; notebooks ship as directly runnable code that PI edits PARAMS in, not templates to copy first.
- The PR plan was rewritten from 10 PRs (one per stage + four infrastructure PRs) down to about 5 PRs (`pyproject` + envs / `read_with_manifest` / `sweep` + scorers + `load_markers` / a batch of stage notebooks / minimal gastric ontology + initial markers). See `_plan.md`.

## Reviewer Cheat Sheet

`code-reviewer` agent uses this table to flag PRs that drift back toward the over-engineered patterns this ADR rolled back. Each row maps a pattern in PR diff to the ADR it violates and the recommended alternative. Reviewer scans this table first; only after no row matches does the agent fall back to general-purpose code review judgement.

| Pattern in PR diff | Violates | Recommended alternative |
|---|---|---|
| New file in `src/scrna_integration/` other than `io.py` / `sweep.py` / `markers.py` / `scorers.py` / `__init__.py` | ADR-0001 + ADR-0004 | Move logic into the relevant stage notebook, or extend an existing module if it's truly within scope. New module = new ADR. |
| `def run_<scanpy_fn>(adata, ...)` wrapping a single `sc.*` call | ADR-0001 | Call `sc.*` directly in the notebook. |
| `@register_panel` / `@register_method` / `@register_scorer` decorator | ADR-0003 | Pass the function directly as an argument (e.g. `custom_panels=[fn1, fn2]`). |
| New parser layer for `manifest` / `config` (re-implementing what `read_with_manifest` already accepts) | ADR-0001 | Caller reads the manifest dict directly with `yaml.safe_load`. |
| `adata.uns["scrna_integration"][...]` reserved namespace write | ADR-0004 | PI's free naming: `adata.uns["my_run_v1"] = {...}`. |
| New `validate_*(adata)` runtime function | ADR-0004 | Code-review enforcement; or a notebook cell with explicit asserts. |
| New framework function added without an ADR | "Three functions" bar (Architectural Stance) | Either drop the new function or write an ADR justifying the bar being met. |
| `si.markers.*` / `si.lineage.*` / `si.tracking.*` / `si.report.*` / `si.qc.*` / `si.disease.*` namespace appearing in imports | ADR-0001 + ADR-0003 + ADR-0004 | These were explicitly removed during the deletion log. Use `from scrna_integration import read_with_manifest, sweep, load_markers` only. |
| YAML / JSON config string referenced by name → looked up in a global dispatch table | ADR-0003 | Pass the actual callable / object directly. |
| New abstraction whose only justification is "extensibility" / "future-proofing" | ADR-0003 | Reject. The bar requires *observed* maintenance pain, not anticipated. |
| `_template_*.ipynb` filename prefix | "Notebooks (Directly Runnable, Not Templates)" section + PI third-round feedback | Notebooks are directly runnable; rename to `stage*.ipynb`. |
| `adata.copy()` without justification in code or PR description | Memory Discipline §2 | Use `inplace=True` (scanpy default); only copy when a separate object is genuinely needed. |
| `adata.write_h5ad(...)` without `compression="lzf"` | Memory Discipline §4 | Always pass `compression="lzf"`. |
| `adata.X` ending up dense after pipeline operations | Memory Discipline §1 | Investigate which scanpy step densified; `adata.X = sp.csr_matrix(adata.X)` to restore. |
| `float64` dtype on `adata.X` post-normalize or on `obsm/X_*` matrices | Memory Discipline §5 | Cast to `float32` after normalize / scvi / harmony output. |
| Stage notebook missing the final `del adata; gc.collect()` cell | Memory Discipline §3 | Add the standard final cell. |
| Stage notebook missing the one-line `assert sp.issparse(adata.X) and adata.X.dtype == np.float32` self-check before `write_h5ad` | Memory Discipline §"In-notebook self-check" | Add the assertion in the cell immediately before `write_h5ad`. |
| Stage notebook drifts from its cell-sequence spec in CONTEXT.md (cells removed / reordered / renamed without justification) | Notebooks §"Cell-sequence specs are binding" | Restore missing/reordered cells, or document the deviation in the PR description. |

Severity follows `code-reviewer`'s standard verdict scale (block / request-changes / approve). Architectural violations (ADR-0001/0003/0004) typically `request-changes` or `block`. Memory-discipline violations typically `request-changes` or appear as advisory comments in an otherwise `approve` verdict.
- Future grilling rounds and PR reviews will treat any new `si.*` namespace as a red flag requiring justification under ADR-0001 / ADR-0003 / this ADR.
