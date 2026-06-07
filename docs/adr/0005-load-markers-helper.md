# `load_markers` as a legitimate third framework function

The framework's `__init__.py` re-exports a third function `load_markers(csv_path, roles=("canonical", "optional"))` for reading the per-tissue marker libraries under `references/markers/`. This is a deliberate exception to ADR-0001 / 0003 / 0004's "framework should not wrap pandas / yaml" stance, justified because PI's prior real-world scRNA-seq workflows confirm the marker-load + role-filter + groupby boilerplate (3–4 lines) recurs across stage 6 annotation, per-cluster profiling, and downstream gene-set scoring notebooks. Centralising the `role` semantics (`canonical` / `optional` / `negative`) prevents the realistic failure mode of a notebook silently including `negative` markers in a gene-set score.

## Considered Options

- **Stay at two functions, accept boilerplate**: rejected. The "≥3 notebooks + naive failure observed" bar from `Architectural Stance` was met by PI's lived experience rather than waiting for it to recur in this project's PR work. The role semantics carry real misuse risk that a 3-line snippet is unlikely to encode safely on every copy-paste.
- **Document the boilerplate in CONTEXT.md but not add a function**: rejected. Documenting "the right way to write this 4-line snippet" still leaves four lines of substantive logic distributed across notebooks; one of those copies will eventually drift.

## Consequences

- The framework's public surface is now three functions: `read_with_manifest`, `sweep`, `load_markers`. CONTEXT.md is updated to reflect this.
- `Architectural Stance > "Two functions" is the current reality` is renamed to "Three functions" and explicitly references this ADR as the precedent for crossing the bar.
- The function lives at `src/scrna_integration/markers.py` and is re-exported from `scrna_integration/__init__.py`. No `markers` sub-namespace is exposed.
- The `roles=None` mode returns the full 3-layer dict, preserving access to all role information for callers that need it (e.g. building a reverse-validation report against `negative` markers).
- Future ADRs proposing a fourth function should still meet the bar described in `Architectural Stance`. This ADR does not relax the bar; it documents one case where it was met.
