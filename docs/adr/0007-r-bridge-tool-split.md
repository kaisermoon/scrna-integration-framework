---
status: accepted (supersedes ADR-0002)
---

# R bridge split by tool: pure-Python first, rpy2 for light in-notebook, subprocess Rscript for heavy R

R access is chosen per tool rather than via one canonical mechanism. The third-round grilling found that the project's real R-using code does not match ADR-0002's "rpy2 for everything" stance: the students' downstream analyses call R via subprocess `Rscript` + temp files (Monocle3, UCell), and several analyses ADR-0002 assumed were R-only are actually pure Python.

## The split

| Scenario | Bridge | Why |
|---|---|---|
| InferCNV, CytoTRACE | **Pure Python** — `infercnvpy` (`cnv.tl.infercnv` / `cnv.tl.cnv_score`), `cellrank` `CytoTRACEKernel`. No R at all. | Mature Python packages exist and are validated in `student-code`. ADR-0002 and an earlier SPEC wrongly tagged InferCNV as rpy2. |
| SoupX (stage 2 ambient correction) | **rpy2 `%%R`** in-notebook | Per-sample matrices are small, conversion is stable, and it keeps QC visually coherent inside one notebook. Validated in `legacy-GCPL/02_quality_control.ipynb`. |
| Monocle3, UCell, DESeq2, hdWGCNA, InferCNV-R variants, other heavy / R-only tools | **subprocess `Rscript` + temp files** (write `.mtx`/`.csv` → `Rscript --vanilla` → read results back) | rpy2 + anndata2ri is brittle when converting heavy objects (Seurat/SCE large matrices) and breaks on numpy/scipy/anndata upgrades. Process isolation is more robust, the `.R` script is independently debuggable, and this is the pattern `student-code` actually uses. |

## Considered Options

- **rpy2 for all R access (ADR-0002)**: rejected. It treats the brittle path as canonical and the robust path as an escape hatch — backwards relative to how the students' working code is written. rpy2's fragility under dependency upgrades is real and already cost the students enough to make them choose subprocess.
- **subprocess Rscript for all R access**: rejected. SoupX's per-sample small-matrix exchange is genuinely cleaner with rpy2 `%%R` inside the QC notebook, and GCPL already runs it that way. Forcing SoupX through temp files would add I/O ceremony for no gain.

## Consequences

- The `notebooks/` R idioms differ by tool: SoupX uses `%load_ext rpy2.ipython` + `%%R -i ... -o ...`; Monocle3/UCell/DESeq2 use a Python cell that writes inputs, calls `subprocess.run(["Rscript", ...])`, and reads outputs back. SPEC's R Bridge section documents both, with subprocess as the primary path for heavy tools.
- Notebook-internal coherence is partially traded away for the heavy-R tools (intermediate temp files, a separate `.R` script, managing the `Rscript` path). The notebook still reads R-produced figures back and displays them, preserving the "everything visible in the notebook" goal at the result level. Accepted by PI.
- `environment-r.yml` still ships the R packages; whether a tool is reached via rpy2 or subprocess does not change which R packages are installed.
- stage7/cnv.ipynb is pure-Python `infercnvpy`, not rpy2. SPEC's stage 7 module table and cnv cell-sequence are corrected accordingly.
- The R env stays structured by stage so users skip R packages they don't need (e.g. a QC-only user installs SoupX but not Monocle3).
