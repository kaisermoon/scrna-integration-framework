# Thin framework over native scanpy

The framework calls scanpy functions directly in notebooks rather than wrapping them in framework classes. PI explicitly rejected the wrapping approach: "scanpy functions cannot be exhaustively enumerated, and wrapping harms readability since everyone already knows scanpy and anndata". The framework adds a thin layer only at four gaps scanpy does not cover (multi-source IO + obs schema, sweep harness, stage reports, run tracking namespace), keeping scanpy as the canonical API surface for ~80% of pipeline operations.

## Considered Options

- **Wrapper classes** (`QCFilter(min_genes=200).fit_transform(adata)`): would unify the interface, integrate naturally with a sweep harness, and provide IDE autocomplete on parameters. Rejected because it forces users to learn a parallel vocabulary alongside scanpy, breaks transferability of scanpy tutorials, and the wrapper layer would need to track the scanpy API surface forever.
- **sklearn-style Pipeline**: rejected for the same reason plus a fundamental fit with scRNA-seq workflow — PI explicitly wants to jump back into earlier stages mid-analysis, which a fixed-pipeline abstraction resists.

## Consequences

- The sweep harness must accept arbitrary callables (including scanpy functions, scvi-tools methods, and user-defined functions) rather than framework-internal step types. Implementation has to be type-flexible.
- Run metadata convention lives by social contract: there is no class-based enforcement that every framework-relevant call writes to `adata.uns["scrna_integration"]`. Code review and stage report assertions enforce it instead.
- Future Stage 2/3 graduation (software paper / PyPI) will not need an API rewrite — the four submodules are already the publishable surface.
