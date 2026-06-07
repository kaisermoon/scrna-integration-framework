"""Sweep harness: iterate over parameter grid, score each combo, write report."""

import itertools
import os
from collections.abc import Callable
from typing import Any

import anndata
import pandas as pd


def sweep(
    fn: Callable,
    adata: anndata.AnnData,
    candidates: dict[str, list[Any]],
    scorer: Callable,
    output_dir: str,
) -> pd.DataFrame:
    """Run *fn* over the Cartesian product of *candidates*, score each, report.

    Args:
        fn: Any callable ``fn(adata, **params)`` — e.g. scanpy or user function.
        adata: Base AnnData; a fresh copy is made for each combination.
        candidates: ``{"param": [v1, v2], ...}`` — the grid to sweep.
        scorer: ``scorer(adata_after, adata_before, params) -> dict[str, float]``.
        output_dir: Directory for the sweep report.

    Returns:
        DataFrame with one row per combo; columns = param names + scorer keys.
    """
    param_names = list(candidates.keys())
    param_value_lists = list(candidates.values())
    rows: list[dict[str, Any]] = []

    os.makedirs(output_dir, exist_ok=True)

    for combo in itertools.product(*param_value_lists):
        params = dict(zip(param_names, combo, strict=False))
        adata_copy = adata.copy()
        fn(adata_copy, **params)
        metrics = scorer(adata_copy, adata, params)
        rows.append({**params, **metrics})

    df = pd.DataFrame(rows)
    _write_report(output_dir, df)
    return df


def _write_report(output_dir: str, df: pd.DataFrame) -> None:
    """Write a minimal markdown report with the sweep results table."""
    lines = ["# Sweep Report\n", f"**{len(df)} combination(s)** evaluated.\n"]

    if df.empty:
        lines.append("_No results._\n")
        with open(os.path.join(output_dir, "sweep_report.md"), "w") as f:
            f.write("\n".join(lines))
        return

    # Separate numeric metrics from string-artifact columns for the table
    metric_cols = [
        c
        for c in df.columns
        if df[c].dtype.kind in ("f", "i") and pd.api.types.is_numeric_dtype(df[c])
    ]
    param_cols = [c for c in df.columns if c not in metric_cols]

    # Build markdown table manually (no tabulate dependency)
    all_cols = param_cols + metric_cols
    header = "| " + " | ".join(all_cols) + " |"
    sep = "|" + "|".join(" --- " for _ in all_cols) + "|"
    lines.extend([header, sep])

    for _, row in df[all_cols].iterrows():
        vals = []
        for col in all_cols:
            v = row[col]
            if isinstance(v, float):
                vals.append(f"{v:.4f}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")

    with open(os.path.join(output_dir, "sweep_report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
