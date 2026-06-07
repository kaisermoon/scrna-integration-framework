"""Marker-library loader for ``references/markers/{tissue}_{purpose}.csv``.

Justified by ADR-0005: the load+filter+groupby boilerplate recurs across
stage-6 annotation, per-cluster profiling, and gene-set scoring notebooks.
Centralising the *role* semantics prevents misuse.
"""

from __future__ import annotations

import pandas as pd


def load_markers(
    csv_path: str,
    roles: tuple[str, ...] | None = ("canonical", "optional"),
) -> dict[str, list[str]] | dict[str, dict[str, list[str]]]:
    """Load a marker CSV and return ``cell_type -> markers``, filtered by role.

    Args:
        csv_path: Any filesystem path (relative or absolute).
        roles: Role filter tuple. Default ``("canonical", "optional")``
            returns the common gene-set-scoring / dotplot case.
            ``roles=("negative",)`` returns only negative markers.
            ``roles=None`` returns the full 3-layer dict
            ``{cell_type: {canonical: [...], optional: [...], negative: [...]}}``.

    Returns:
        Dict keyed by cell type — flat list (when *roles* is a tuple) or
        nested dict (when *roles* is ``None``).

    CSV schema: ``tissue, cell_type, marker, role, reference, notes``
    """
    df = pd.read_csv(csv_path)

    if roles is None:
        result: dict[str, dict[str, list[str]]] = {}
        for role in ("canonical", "optional", "negative"):
            role_df = df[df["role"] == role]
            for ct, group in role_df.groupby("cell_type"):
                result.setdefault(ct, {}).setdefault(role, [])
                result[ct][role] = group["marker"].tolist()
        # Ensure every cell_type has all three role keys
        for ct in result:
            for role in ("canonical", "optional", "negative"):
                result[ct].setdefault(role, [])
        return result

    # Flat mode: filter by role, return {cell_type: [markers]}
    filtered = df[df["role"].isin(roles)]
    return {
        ct: group["marker"].tolist() for ct, group in filtered.groupby("cell_type")
    }
