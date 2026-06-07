---
status: accepted
---

# Absorb student-code downstream techniques by re-implementing to project conventions, not by copying

The third-round grilling reviewed `references/student-code/` and found a broad set of downstream techniques that GCPL lacks and that the project wants: transcriptome entropy, CytoTRACE, multi-metric root-cluster identification, Monocle3 trajectory, UCell gene-set scoring with transition/mixed detection, scCODA compositional analysis, effect-size statistics, gene-existence checking before marker use. PI's decision: **all of these are absorbed into the roadmap**, but each is **re-implemented to this project's conventions**, not lifted verbatim from the student scripts.

The student code is a correctness reference (the algorithm works, the library calls are right) and a source of specific good ideas — not a code donor. Copying it wholesale would import its anti-patterns: hard-coded Windows paths, `!pip install` cells, 800-line monolithic scripts, and the copy-paste-per-cluster template family (the CD4/CD8 deep-analysis scripts are 17 near-identical files — exactly the duplication this project rejects).

## What "re-implement to our conventions" means concretely

A PR-5+ module that absorbs a student technique must:

- Follow the stage-notebook structure (PARAMS cell, scanpy-native body, `adata.uns` run-metadata, memory self-check, `del adata; gc.collect()` ending).
- Use the project's R-bridge split (ADR-0007): pure-Python where a Python package exists (CytoTRACE via cellrank, entropy in numpy), subprocess `Rscript` for heavy R (Monocle3, UCell, scCODA's R parts if any).
- Read markers via `load_markers`; check gene existence against `var_names` before scoring.
- Stay inside the three-function framework boundary (ADR-0001/0003/0004): no new framework functions, no plugin registry, no per-cluster script templating. A reusable helper that genuinely recurs goes through the ADR-0004 escape-hatch bar, not in by default.
- Reference the source explicitly in the PR description: "re-implements the entropy/root-finding logic from `student-code/workflow_for_pseudotime/4.3_*.py` and `4.4_*.py`, rewritten to stage-7 notebook conventions."

## Considered Options

- **Copy the working student scripts as-is, clean up later**: rejected. "Clean up later" rarely happens; the anti-patterns become load-bearing. The 项目构思 explicitly wants simple, well-annotated, scanpy-native code for PI + non-CS students — the student monoliths are the opposite.
- **Only reference, don't absorb (leave downstream to ad-hoc per-project scripts)**: rejected. PI wants these techniques as standing, reusable stage-7 modules in the framework's notebook set, not re-derived each project.

## Consequences

- SPEC's stage-7 module table gains a "student-code reference" column pointing each module at the file(s) to re-implement from.
- The PR plan phases stage-7 expansion one module per PR (pseudotime, abundance, pathway, GRN, cell-communication, gene-modules), each absorbing its student techniques. No single PR carries the whole downstream.
- The CD4/CD8 deep-analysis template family is **reference-only** and is never templatized into the framework; subset deep-dives use stage 6.5 + the relevant stage-7 module notebooks.
- gene-existence checking and UCell transition/mixed detection are cross-cutting idioms documented where used (marker convention; gene-set scoring cell), not separate framework functions.
