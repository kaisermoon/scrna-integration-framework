"""Tests for sweep() harness — synthetic AnnData, no real data files."""

import os
import tempfile

import anndata
import numpy as np
import pandas as pd

from scrna_integration.sweep import sweep


def _dummy_adata(n_cells: int = 100, n_genes: int = 50) -> anndata.AnnData:
    rng = np.random.default_rng(42)
    X = rng.poisson(1, size=(n_cells, n_genes)).astype(np.float32)  # noqa: N806
    adata = anndata.AnnData(X)
    adata.obs["n_genes"] = rng.integers(200, 5000, size=n_cells).astype(np.float64)
    return adata


def _filter_fn(adata: anndata.AnnData, threshold: float = 500, **kwargs) -> anndata.AnnData:
    """Dummy: filter cells by n_genes threshold. **kwargs absorbs unused sweep params."""
    adata._inplace_subset_obs(adata.obs["n_genes"] > threshold)
    return adata


def _count_scorer(adata_after, adata_before, params):
    return {"cells_remaining": adata_after.n_obs, "retention": adata_after.n_obs / max(adata_before.n_obs, 1)}


class TestSweep:
    def test_cartesian_product_combinations(self):
        """sweep generates the correct number of Cartesian product combos."""
        adata = _dummy_adata(100)
        candidates = {"threshold": [200, 500, 1000], "extra": ["a", "b"]}
        df = sweep(_filter_fn, adata, candidates, _count_scorer, tempfile.mkdtemp())
        assert len(df) == 6

    def test_returns_dataframe_with_correct_columns(self):
        adata = _dummy_adata(100)
        candidates = {"threshold": [500]}
        df = sweep(_filter_fn, adata, candidates, _count_scorer, tempfile.mkdtemp())
        assert isinstance(df, pd.DataFrame)
        assert "threshold" in df.columns
        assert "cells_remaining" in df.columns
        assert "retention" in df.columns

    def test_does_not_mutate_original_adata(self):
        """Original adata is untouched after sweep."""
        adata = _dummy_adata(100)
        original_n_obs = adata.n_obs
        sweep(_filter_fn, adata, {"threshold": [500]}, _count_scorer, tempfile.mkdtemp())
        assert adata.n_obs == original_n_obs

    def test_report_generated(self):
        adata = _dummy_adata(100)
        out = tempfile.mkdtemp()
        sweep(_filter_fn, adata, {"threshold": [500]}, _count_scorer, out)
        rp = os.path.join(out, "sweep_report.md")
        assert os.path.exists(rp)
        content = open(rp).read()
        assert "threshold" in content or "Sweep" in content

    def test_wraps_arbitrary_callable(self):
        """sweep wraps any callable, not just scanpy functions."""

        def custom_fn(adata, multiplier=2):
            adata.X = adata.X * multiplier
            return adata

        def custom_scorer(adata_after, adata_before, params):
            return {"mean_x": float(np.mean(adata_after.X))}

        adata = _dummy_adata(50, 20)
        df = sweep(custom_fn, adata, {"multiplier": [2, 3]}, custom_scorer, tempfile.mkdtemp())
        assert len(df) == 2
        assert "mean_x" in df.columns

    def test_single_candidate_value(self):
        adata = _dummy_adata(100)
        df = sweep(_filter_fn, adata, {"threshold": [500]}, _count_scorer, tempfile.mkdtemp())
        assert len(df) == 1

    def test_empty_candidates(self):
        """Empty candidates dict runs once with no params."""
        adata = _dummy_adata(100)
        df = sweep(_filter_fn, adata, {}, _count_scorer, tempfile.mkdtemp())
        assert len(df) == 1
