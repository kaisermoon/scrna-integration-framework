---
status: superseded by ADR-0007
---

# R interoperability via rpy2 + anndata2ri

> **Superseded by [ADR-0007](./0007-r-bridge-tool-split.md).** The "rpy2 for all R access, subprocess as rare escape hatch" stance below was reversed after the third-round grilling examined `references/student-code/`: the students' real downstream R calls (Monocle3, UCell) are all subprocess `Rscript`, not rpy2. ADR-0007 splits R access by tool. The `_r_bridge` module referenced here was also already removed by the framework thinning (ADR-0004) — there is no centralised rpy2 module. This ADR is retained for history.

The framework uses `rpy2` + `anndata2ri` for all R access — IO (Seurat/SingleCellExperiment readers), QC (SoupX), and downstream stages (Monocle3, InferCNV, SCENIC, scran, DESeq2). The alternative — wrapping each R algorithm in an independent `Rscript` subprocess — was rejected because rpy2 is already a hard dependency for QC stage 2 (SoupX in legacy GCPL code), so a subprocess path would add a second parallel R-call mechanism rather than removing the rpy2 dependency.

## Considered Options

- **Subprocess `Rscript` per algorithm**: rejected per above. Originally proposed to avoid rpy2's known fragility under numpy/scipy/anndata upgrades, but the trade-off is upside-down once rpy2 is unavoidable elsewhere.
- **rpy2 only, no anndata2ri**: rejected — manual `SingleCellExperiment ↔ AnnData` conversion is too error-prone given the matrix orientation flips and obs/var slot remapping involved.

## Consequences

- Subprocess `_r_bridge.run_r_script(...)` is preserved as an escape hatch for the case where rpy2 conversion breaks after dependency upgrades. Use is rare and documented per incident.
- The framework ships two conda environment files (`environment.yml` Python, `environment-r.yml` R). The R environment is structured by stage so users can install only the R packages they currently need.
- All rpy2 boilerplate is centralised in `scrna_integration.io._r_bridge`. Stages must not roll their own `from rpy2.robjects import ...` blocks.
- Future Stage 2/3 graduation (PyPI) will need to document the dual-environment installation prominently — this is friction we accept for now.
