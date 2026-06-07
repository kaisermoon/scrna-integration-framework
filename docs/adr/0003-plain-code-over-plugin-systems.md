# Plain code over plugin systems

When the framework needs to accept user-provided extensions — custom report panels, custom scorers, custom QC steps, custom annotation methods — the default solution is a plain function (or list of plain functions) passed as an argument. Registries, decorators, plugin discovery, dynamic imports, and similar "framework magic" are rejected unless a plain-function approach has been demonstrably tried and failed.

This rule generalises ADR-0001 (thin framework over scanpy). The same reasoning applies: PI, students, and reviewing agents all read framework code routinely; obscuring the call graph behind decorator-based registration costs more in readability than it saves in API ergonomics. The framework's audience is small and known, so user-extensibility patterns designed for open-source libraries (scikit-learn pluggable estimators, pytest fixture discovery, etc.) are over-built for this context.

## What this looks like in practice

- **Report panels**: `custom_panels: list[Callable]` rather than `@register_panel`.
- **Scorers**: `scorers: dict[str, Callable]` rather than scorer plugins.
- **QC steps**: explicit step list in stage entry function rather than discoverable QC plugin namespace.
- **New methods**: covered by Sweep harness's "any callable" contract — no `@register_method` decorator.

## Considered Options

- **Registry / decorator pattern**: rejected as the default. Allowed only when call sites genuinely cannot reference the implementation by name (e.g. plugins loaded from outside this repo, which the framework does not need).
- **Configuration-driven dispatch** (YAML lists method names looked up via string): rejected. Stringly-typed lookups break IDE navigation, type checking, and refactor-safety.

## Consequences

- Code review (`code-reviewer` agent) treats registries and decorators as red flags requiring explicit justification.
- New PRs that add framework-internal abstractions claiming "extensibility" will be challenged to demonstrate why a plain function passed as argument cannot serve.
- Documentation pages on "how to add a custom X" become trivial — they show a normal function definition and a list-passing call site, not a registration ceremony.
- The framework's `__all__` / public API stays small. There is no plugin enumeration logic to maintain.
