---
status: accepted (supersedes ADR-0002; Soupx subprocess revision 2026-06-09)
---

# R bridge split by tool: pure-Python first, subprocess Rscript for all R (rpy2 abandoned)

R access is chosen per tool rather than via one canonical mechanism. The third-round grilling found that the project's real R-using code does not match ADR-0002's "rpy2 for everything" stance: the students' downstream analyses call R via subprocess `Rscript` + temp files (Monocle3, UCell), and several analyses ADR-0002 assumed were R-only are actually pure Python.

**2026-06-09 修订**: SoupX 从 rpy2 `%%R` 改判为 subprocess Rscript（见下方"SoupX 改判"节）。至此，项目的 R 桥接策略统一为 subprocess Rscript 全覆盖，不再保留 rpy2 路径。

## The split

| Scenario | Bridge | Why |
|---|---|---|
| InferCNV, CytoTRACE | **Pure Python** — `infercnvpy` (`cnv.tl.infercnv` / `cnv.tl.cnv_score`), `cellrank` `CytoTRACEKernel`. No R at all. | Mature Python packages exist and are validated in `student-code`. ADR-0002 and an earlier SPEC wrongly tagged InferCNV as rpy2. |
| SoupX (stage 2 ambient correction), Monocle3, UCell, DESeq2, hdWGCNA, and all other R tools | **subprocess `Rscript` + temp files** (write `.mtx`/`.csv` → `Rscript --vanilla` → read results back) | Process isolation is more robust, the `.R` script is independently debuggable, and this is the pattern `student-code` actually uses. rpy2 + anndata2ri was attempted for SoupX but proved infeasible (see "SoupX reclassification" below). |

## Considered Options

- **rpy2 for all R access (ADR-0002)**: rejected. It treats the brittle path as canonical and the robust path as an escape hatch — backwards relative to how the students' working code is written. rpy2's fragility under dependency upgrades is real and already cost the students enough to make them choose subprocess.
- **subprocess Rscript for all R access**: **accepted (2026-06-09 revision)**. Originally rejected for SoupX due to perceived overhead for small-matrix exchange; the original ADR kept SoupX on rpy2 for "QC coherence inside one notebook." This proved infeasible in practice (see "SoupX reclassification" below). With SoupX moved to subprocess, the split collapses to a single strategy: all R tools go through subprocess Rscript.
- **rpy2 for light tools, subprocess for heavy (original ADR-0007)**: superseded by the 2026-06-09 revision. The "light" path (rpy2 for SoupX) could not be made to work.

### SoupX reclassification (2026-06-09, PI decision B)

**Original classification**: rpy2 `%%R` in-notebook. Rationale at the time: per-sample matrices are small, conversion is stable, keeps QC visually coherent inside one notebook.

**Why reclassified**:
1. **rpy2 + anndata2ri 桥接物理不可行**。conda R 4.4.3 环境下 `R_getVar` 符号缺失（`ImportError: cannot import name 'R_getVar'`）。系统 R 4.2.3 同样缺失。该符号来自已删除的 R framework，pip 安装的 rpy2 链接了不存在的库路径，无法通过重装 rpy2 或 R 版本解决。
2. **subprocess 工艺隔离更稳**。与 DESeq2 / Monocle3 / hdWGCNA 等重型 R 工具统一模式——Python 侧导出矩阵到临时文件、subprocess 调 Rscript 执行独立 `.R` 脚本、读回结果。R 脚本可脱离 notebook 独立调试（`Rscript --vanilla scripts/soupx_run.R <args>`）。
3. **notebook 可读性不降反升**。原 rpy2 `%%R` cell 混合 Python 和 R 语法，非 CS 学生阅读困难。新模式下 Python cell 自包含（导出→调用→读回），R 逻辑集中在 `scripts/soupx_run.R`，两边各自独立可读。

**代价**：SoupX 临时文件落入 `results/_soupx_tmp/`（gitignored），与 DESeq2 的 `results/_deseq2_tmp/` 同模式。PI 接受此代价。

## Consequences

- **All R tools use the same subprocess pattern**. SoupX (`scripts/soupx_run.R`), Monocle3, UCell, DESeq2 (`scripts/deseq2_contrast.R`), hdWGCNA — all follow: Python writes inputs (`.mtx`/`.csv`) → `subprocess.run(["Rscript", "--vanilla", ...])` → reads outputs back. No rpy2 dependency anywhere in the pipeline (except `read_with_manifest` handling `format: "rds"`).
- **`.R` scripts are independently debuggable**. Each lives in `scripts/` and can be run standalone: `Rscript --vanilla scripts/soupx_run.R <args>`. This makes R-side debugging faster and decoupled from the Python notebook.
- **Notebook-internal coherence at the result level**. Intermediate temp files, a separate `.R` script, and managing the `Rscript` path are accepted tradeoffs. The notebook still displays R-produced results inline, preserving "everything visible in the notebook" at the result level.
- `environment-r.yml` still ships the R packages; the bridge mechanism change does not affect which R packages are installed.
- The R env stays structured by stage so users skip R packages they don't need (e.g. a QC-only user installs SoupX but not Monocle3).
