"""Per-leaf static JSON/AST tests for P1-c-template doublet three-state logic.

Tests validate the 01_template_10x.ipynb notebook without executing cells
(no manifest/real data needed).  Follows the test_pr1c1_stage04.py pattern
of json.loads + marker-anchored cell extraction.
"""

import ast
import json
import re
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks/01_per_dataset/01_template_10x.ipynb"


def _source(cell: dict) -> str:
    value = cell.get("source", "")
    return "".join(value) if isinstance(value, list) else value


def _cell(marker: str) -> str:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return next(
        _source(cell)
        for cell in notebook["cells"]
        if cell["cell_type"] == "code" and marker in _source(cell)
    )


def _all_code_sources() -> list[str]:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return [
        "".join(c["source"]) if isinstance(c["source"], list) else c["source"]
        for c in notebook["cells"]
        if c["cell_type"] == "code"
    ]


class TestNotebookSyntax:
    """Validate the notebook can be parsed and all code cells are syntactically valid."""

    def test_notebook_json_parsable(self) -> None:
        notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        assert isinstance(notebook, dict)
        assert "cells" in notebook

    def test_all_code_cells_parse_individually(self) -> None:
        nb = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        for i, cell in enumerate(nb["cells"]):
            if cell["cell_type"] != "code":
                continue
            src = _source(cell)
            cid = cell.get("id", "?")
            try:
                ast.parse(src, filename=f"cell_{i}_{cid}")
            except SyntaxError as e:
                pytest.fail(f"cell[{i}] id={cid}: {e}")


class TestParamsCell:
    """Validate PARAMS cell contains the three new P1-c parameters."""

    def test_new_doublet_params_present(self) -> None:
        params = _cell("# === PARAMS（按数据集修改此 cell）===")
        for var in ["DOUBLET_UNCERTAIN_MARGIN", "DOUBLET_MIN_CELLS", "DOUBLET_MAX_RATE_WARN"]:
            assert var in params, f"PARAMS cell missing {var}"

    def test_default_values_unchanged(self) -> None:
        params = _cell("# === PARAMS（按数据集修改此 cell）===")
        # Extract assignment values via simple regex (not full eval)
        m = re.search(r"DOUBLET_UNCERTAIN_MARGIN\s*=\s*([0-9.]+)", params)
        assert m is not None, "DOUBLET_UNCERTAIN_MARGIN assignment not found"
        assert float(m.group(1)) == 0.2, f"Expected 0.2, got {m.group(1)}"

        m = re.search(r"DOUBLET_MIN_CELLS\s*=\s*([0-9]+)", params)
        assert m is not None, "DOUBLET_MIN_CELLS assignment not found"
        assert int(m.group(1)) == 50, f"Expected 50, got {m.group(1)}"

        m = re.search(r"DOUBLET_MAX_RATE_WARN\s*=\s*([0-9.]+)", params)
        assert m is not None, "DOUBLET_MAX_RATE_WARN assignment not found"
        assert float(m.group(1)) == 0.30, f"Expected 0.30, got {m.group(1)}"

    def test_original_params_retained(self) -> None:
        params = _cell("# === PARAMS（按数据集修改此 cell）===")
        for var in ["EXPECTED_DOUBLET_RATE", "DOUBLET_SCORE_THRESHOLD", "RANDOM_SEED"]:
            assert var in params, f"PARAMS cell missing original param {var}"


class TestScrubletCell:
    """Validate the rewritten scrublet cell contains three-state logic."""

    def test_three_state_literals_present(self) -> None:
        scrub = _cell("三态：singlet / uncertain / doublet")
        for state in ['"singlet"', '"uncertain"', '"doublet"']:
            assert state in scrub, f"Scrublet cell missing literal {state}"

    def test_new_obs_columns_assigned(self) -> None:
        scrub = _cell("三态：singlet / uncertain / doublet")
        for col in ["doublet_prediction", "doublet_include", "doublet_threshold"]:
            assert f'obs["{col}"]' in scrub or f"obs[{repr(col)}]" in scrub, (
                f"Scrublet cell missing obs column assignment for {col}"
            )

    def test_doublet_detection_uns_written(self) -> None:
        scrub = _cell("三态：singlet / uncertain / doublet")
        assert 'adata.uns["doublet_detection"]' in scrub, (
            "Scrublet cell missing uns doublet_detection write"
        )

    def test_no_old_predicted_doublet_assignment(self) -> None:
        scrub = _cell("三态：singlet / uncertain / doublet")
        pattern = re.compile(r'obs\["predicted_doublet"\]\s*=')
        assert not pattern.search(scrub), (
            "Scrublet cell contains obsolete predicted_doublet boolean assignment"
        )

    def test_doublet_needs_review_module_var_set(self) -> None:
        scrub = _cell("三态：singlet / uncertain / doublet")
        assert "DOUBLET_NEEDS_REVIEW" in scrub, (
            "Scrublet cell missing DOUBLET_NEEDS_REVIEW module-level variable"
        )

    def test_per_sample_copy_and_del_gc(self) -> None:
        scrub = _cell("三态：singlet / uncertain / doublet")
        assert "sub = adata[" in scrub, "Missing per-sample adata copy"
        assert "del sub" in scrub, "Missing del sub memory cleanup"
        assert "gc.collect()" in scrub, "Missing gc.collect()"


class TestNoSubsetIn01:
    """Validate that 01 never physically subsets by doublet columns."""

    def test_no_doublet_based_subset(self) -> None:
        all_sources = _all_code_sources()
        all_code = "\n".join(all_sources)
        pattern = re.compile(
            r"adata\s*=\s*adata\s*\[.*doublet_(include|prediction)"
        )
        assert not pattern.search(all_code), (
            "01 notebook contains doublet-based subset "
            "(exclusion must happen in 02)"
        )

    def test_mad_filter_cell_no_doublet_subset(self) -> None:
        mad_cell = _cell("# === 应用 QC 阈值 ===")
        pattern = re.compile(r"doublet_(include|prediction)")
        assert not pattern.search(mad_cell), (
            "MAD filter cell references doublet columns "
            "(should only filter by QC metrics)"
        )


class TestAssertionCell:
    """Validate the assertion cell contains three-state contract checks."""

    def test_doublet_prediction_set_assertion(self) -> None:
        assertions = _cell("# === 断言对齐契约（进入 02_merged 前的硬闸门）===")
        assert 'set(adata.obs["doublet_prediction"]' in assertions, (
            "Assertion cell missing doublet_prediction set check"
        )
        assert '"singlet", "uncertain", "doublet"' in assertions, (
            "Assertion cell missing explicit three-state values"
        )

    def test_doublet_include_dtype_check(self) -> None:
        assertions = _cell("# === 断言对齐契约（进入 02_merged 前的硬闸门）===")
        assert 'doublet_include" in adata.obs' in assertions, (
            "Assertion cell missing doublet_include existence check"
        )
        assert 'doublet_include"].dtype == bool' in assertions, (
            "Assertion cell missing doublet_include dtype bool check"
        )

    def test_qc_report_three_state_fields(self) -> None:
        assertions = _cell("# === 断言对齐契约（进入 02_merged 前的硬闸门）===")
        for field in ["n_doublet", "n_uncertain", "n_included", "doublet_needs_review"]:
            assert field in assertions, (
                f"qc_report missing three-state field: {field}"
            )

    def test_qc_report_no_old_predicted_doublet(self) -> None:
        assertions = _cell("# === 断言对齐契约（进入 02_merged 前的硬闸门）===")
        assert "predicted_doublet" not in assertions or '"n_doublet_predicted"' not in assertions, (
            "qc_report still references obsolete predicted_doublet"
        )


class TestDiagnosticCell:
    """Validate the per-sample doublet diagnostic cell exists."""

    def test_diagnostic_cell_exists(self) -> None:
        nb = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        diag_cells = [
            c for c in nb["cells"]
            if c["cell_type"] == "code"
            and 'doublet_detection"]["per_sample_summary"' in _source(c)
        ]
        assert len(diag_cells) >= 1, (
            "Per-sample doublet diagnostic code cell not found"
        )

    def test_diagnostic_include_impact(self) -> None:
        diag = _cell('doublet_detection"]["per_sample_summary"')
        assert "n_incl" in diag and "n_excl" in diag, (
            "Diagnostic cell missing included/excluded count"
        )
        assert "median_n_genes" in diag, (
            "Diagnostic cell missing QC comparison"
        )


class TestThreeStateBoundaryLogic:
    """Isolated unit test for the three-state classification boundary.

    These tests verify the classification logic in isolation by simulating
    the three-state decision with numpy masks (matching the notebook pattern).
    They do NOT require a running scrublet kernel.
    """

    @staticmethod
    def _classify(
        scores: "np.ndarray",
        threshold_high: float,
        uncertain_margin: float,
    ) -> "tuple[np.ndarray, np.ndarray]":
        """Reproduce the notebook's three-state classification logic.

        Returns (predictions, include_mask).
        """
        threshold_low = threshold_high * (1 - uncertain_margin)
        pred = np.full(len(scores), "singlet", dtype=object)
        pred[scores >= threshold_high] = "doublet"
        pred[(scores >= threshold_low) & (scores < threshold_high)] = "uncertain"
        include = pred != "doublet"
        return pred, include

    def test_score_at_high_threshold_is_doublet(self) -> None:
        pred, include = self._classify(
            np.array([0.30]), threshold_high=0.30, uncertain_margin=0.2,
        )
        assert pred[0] == "doublet"
        assert not include[0]

    def test_score_at_low_threshold_is_uncertain(self) -> None:
        pred, include = self._classify(
            np.array([0.24]), threshold_high=0.30, uncertain_margin=0.2,
        )
        # threshold_low = 0.30 * (1 - 0.2) = 0.24
        assert pred[0] == "uncertain"
        # uncertain still included
        assert include[0]

    def test_score_below_low_threshold_is_singlet(self) -> None:
        pred, include = self._classify(
            np.array([0.23]), threshold_high=0.30, uncertain_margin=0.2,
        )
        assert pred[0] == "singlet"
        assert include[0]

    def test_include_mask_equals_not_doublet(self) -> None:
        rng = np.random.default_rng(42)
        scores = rng.uniform(0, 0.5, size=100)
        pred, include = self._classify(scores, threshold_high=0.25, uncertain_margin=0.2)
        expected_include = pred != "doublet"
        assert (include == expected_include).all(), (
            "doublet_include must equal (prediction != 'doublet')"
        )

    def test_zero_uncertain_margin_no_uncertain(self) -> None:
        pred, include = self._classify(
            np.array([0.20, 0.24999, 0.25, 0.26]),
            threshold_high=0.25, uncertain_margin=0.0,
        )
        # threshold_low = threshold_high, no uncertain zone
        assert list(pred) == ["singlet", "singlet", "doublet", "doublet"]

    def test_singlet_uncertain_doublet_are_exhaustive(self) -> None:
        rng = np.random.default_rng(99)
        scores = rng.uniform(0, 0.5, size=200)
        pred, _ = self._classify(scores, threshold_high=0.25, uncertain_margin=0.2)
        unique = set(pred)
        assert unique <= {"singlet", "uncertain", "doublet"}, (
            f"Unexpected prediction values: {unique}"
        )


class TestObsAndUnsSchema:
    """Validate that the field naming matches the P1-c queue specification."""

    def test_obs_column_names_in_code(self) -> None:
        """Check all four required obs columns appear in the notebook code."""
        all_sources = _all_code_sources()
        all_code = "\n".join(all_sources)
        for col in ["doublet_score", "doublet_prediction", "doublet_threshold", "doublet_include"]:
            assert f'obs["{col}"]' in all_code or f"{repr(col)}" in all_code, (
                f"Required obs column {col} not found in notebook code"
            )

    def test_uns_key_names_in_code(self) -> None:
        """Check that doublet_detection uns key and its sub-keys appear."""
        all_sources = _all_code_sources()
        all_code = "\n".join(all_sources)
        assert 'uns["doublet_detection"]' in all_code, (
            "Required uns key doublet_detection not found"
        )
        for sub_key in ["per_sample_summary", "needs_review", "skip", "skip_reason", "params"]:
            assert sub_key in all_code, (
                f"Required uns sub-key {sub_key} not found in notebook code"
            )

    def test_qc_report_fields_in_code(self) -> None:
        assertions = _cell("# === 断言对齐契约（进入 02_merged 前的硬闸门）===")
        for field in ["n_doublet", "n_uncertain", "n_included", "doublet_needs_review"]:
            assert f'"{field}"' in assertions or f"'{field}'" in assertions, (
                f"qc_report missing field: {field}"
            )
