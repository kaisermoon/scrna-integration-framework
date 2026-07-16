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
    """Validate PARAMS cell contains the new doublet three-state parameters."""

    PARAMS_MARKER = "# === PARAMS：本 cell 是唯一修改参数的地方 ==="

    def test_new_doublet_params_present(self) -> None:
        params = _cell(self.PARAMS_MARKER)
        for var in ["DOUBLET_SCORE_HIGH", "DOUBLET_SCORE_LOW",
                     "DOUBLET_UNCERTAIN_INCLUDE", "DOUBLET_RATE_ALERT_HIGH",
                     "DOUBLET_MIN_CELLS"]:
            assert var in params, f"PARAMS cell missing {var}"

    def test_default_values_unchanged(self) -> None:
        params = _cell(self.PARAMS_MARKER)
        # DOUBLET_SCORE_HIGH = None
        m = re.search(r"DOUBLET_SCORE_HIGH\s*=\s*(None|[0-9.]+)", params)
        assert m is not None, "DOUBLET_SCORE_HIGH assignment not found"
        assert m.group(1) == "None", f"Expected None, got {m.group(1)}"

        # DOUBLET_SCORE_LOW = None
        m = re.search(r"DOUBLET_SCORE_LOW\s*=\s*(None|[0-9.]+)", params)
        assert m is not None, "DOUBLET_SCORE_LOW assignment not found"
        assert m.group(1) == "None", f"Expected None, got {m.group(1)}"

        # DOUBLET_UNCERTAIN_INCLUDE = True
        m = re.search(r"DOUBLET_UNCERTAIN_INCLUDE\s*=\s*(True|False)", params)
        assert m is not None, "DOUBLET_UNCERTAIN_INCLUDE assignment not found"
        assert m.group(1) == "True", f"Expected True, got {m.group(1)}"

        # DOUBLET_RATE_ALERT_HIGH = 0.40
        m = re.search(r"DOUBLET_RATE_ALERT_HIGH\s*=\s*([0-9.]+)", params)
        assert m is not None, "DOUBLET_RATE_ALERT_HIGH assignment not found"
        assert float(m.group(1)) == 0.40, f"Expected 0.40, got {m.group(1)}"

        # DOUBLET_MIN_CELLS = 50
        m = re.search(r"DOUBLET_MIN_CELLS\s*=\s*([0-9]+)", params)
        assert m is not None, "DOUBLET_MIN_CELLS assignment not found"
        assert int(m.group(1)) == 50, f"Expected 50, got {m.group(1)}"

    def test_original_params_retained(self) -> None:
        params = _cell(self.PARAMS_MARKER)
        for var in ["EXPECTED_DOUBLET_RATE", "DOUBLET_SCORE_THRESHOLD", "RANDOM_SEED"]:
            assert var in params, f"PARAMS cell missing original param {var}"


class TestScrubletCell:
    """Validate the rewritten scrublet cell contains three-state logic."""

    SCRUBLET_MARKER = "# === 双细胞鉴定：per-sample scrublet（三态定级）==="

    def test_three_state_literals_present(self) -> None:
        scrub = _cell(self.SCRUBLET_MARKER)
        for state in ['"singlet"', '"uncertain"', '"doublet"']:
            assert state in scrub, f"Scrublet cell missing literal {state}"

    def test_new_obs_columns_assigned(self) -> None:
        scrub = _cell(self.SCRUBLET_MARKER)
        for col in ["doublet_score", "doublet_class", "doublet_include"]:
            assert f'obs["{col}"]' in scrub or f"obs[{repr(col)}]" in scrub, (
                f"Scrublet cell missing obs column assignment for {col}"
            )

    def test_doublet_contract_uns_written(self) -> None:
        # doublet_contract is written in the 收尾 cell (22), not in the per-sample
        # scrublet cell itself.  Verify the contract pattern exists somewhere.
        all_code = "\n".join(_all_code_sources())
        assert 'adata.uns["doublet_contract"]' in all_code, (
            "Notebook missing uns doublet_contract write"
        )

    def test_predicted_doublet_derived_from_doublet_class(self) -> None:
        """predicted_doublet is now a backward-compat boolean derived from doublet_class."""
        all_code = "\n".join(_all_code_sources())
        assert 'predicted_doublet' in all_code, (
            "Notebook missing backward-compat predicted_doublet"
        )
        assert 'doublet_class' in all_code, (
            "Notebook missing primary doublet_class column"
        )

    def test_doublet_needs_review_var_set(self) -> None:
        scrub = _cell(self.SCRUBLET_MARKER)
        assert "doublet_needs_review" in scrub, (
            "Scrublet cell missing doublet_needs_review local variable"
        )

    def test_per_sample_loop_structure(self) -> None:
        scrub = _cell(self.SCRUBLET_MARKER)
        assert "for sample_id in sorted(" in scrub, "Missing per-sample loop"
        assert "sub = adata[" in scrub, "Missing per-sample adata copy"
        assert "sample_mask" in scrub, "Missing sample_mask"


class TestNoSubsetIn01:
    """Validate that 01 never physically subsets by doublet columns."""

    def test_no_doublet_based_subset(self) -> None:
        all_sources = _all_code_sources()
        all_code = "\n".join(all_sources)
        pattern = re.compile(
            r"adata\s*=\s*adata\s*\[.*doublet_(include|class|prediction)"
        )
        assert not pattern.search(all_code), (
            "01 notebook contains doublet-based subset "
            "(exclusion must happen in 02)"
        )

    def test_qc_filter_cell_no_doublet_subset(self) -> None:
        qc_cell = _cell("# === QC 过滤：应用阈值 ===")
        pattern = re.compile(r"doublet_(include|class|prediction)")
        assert not pattern.search(qc_cell), (
            "QC filter cell references doublet columns "
            "(should only filter by QC metrics)"
        )


class TestAssertionCell:
    """Validate the checkpoint cell contains doublet postcondition checks."""

    CHECKPOINT_MARKER = "# === Checkpoint：run_contract + 写入 h5ad + schema 校验 ==="

    def test_doublet_class_set_assertion(self) -> None:
        checkpoint = _cell(self.CHECKPOINT_MARKER)
        assert 'set(adata.obs["doublet_class"]' in checkpoint, (
            "Checkpoint missing doublet_class set check"
        )
        assert '"singlet", "uncertain", "doublet"' in checkpoint, (
            "Checkpoint missing explicit three-state values"
        )

    def test_doublet_include_consistency_check(self) -> None:
        checkpoint = _cell(self.CHECKPOINT_MARKER)
        assert '"doublet_include"' in checkpoint, (
            "Checkpoint missing doublet_include reference"
        )
        assert 'doublet_include_consistent' in checkpoint, (
            "Checkpoint missing doublet_include_consistent postcondition"
        )

    def test_qc_report_three_state_fields(self) -> None:
        qc_report_cell = _cell("# === QC 报告摘要 ===")
        for field in ["n_singlet", "n_uncertain", "n_doublet", "n_excluded"]:
            assert field in qc_report_cell, (
                f"qc_report missing field: {field}"
            )

    def test_qc_report_no_old_predicted_doublet(self) -> None:
        qc_report_cell = _cell("# === QC 报告摘要 ===")
        assert "predicted_doublet" not in qc_report_cell, (
            "qc_report still references obsolete predicted_doublet"
        )


class TestDiagnosticCell:
    """Validate the per-sample doublet diagnostic information exists."""

    def test_diagnostic_contract_present(self) -> None:
        nb = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        diag_cells = [
            c for c in nb["cells"]
            if c["cell_type"] == "code"
            and '"per_sample_diagnostics"' in _source(c)
        ]
        assert len(diag_cells) >= 1, (
            "Per-sample doublet diagnostic code (per_sample_diagnostics) not found"
        )

    def test_diagnostic_excluded_count(self) -> None:
        diag = _cell("# === 双细胞三态收尾：include 派生 + doublet_contract 持久化 ===")
        assert "n_excluded" in diag, (
            "Diagnostic cell missing excluded count"
        )
        assert "n_doublet" in diag or "n_singlet" in diag, (
            "Diagnostic cell missing doublet/singlet counts"
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
        """Check all required obs columns appear in the notebook code."""
        all_sources = _all_code_sources()
        all_code = "\n".join(all_sources)
        for col in ["doublet_score", "doublet_class", "doublet_include"]:
            assert f'obs["{col}"]' in all_code or f"{repr(col)}" in all_code, (
                f"Required obs column {col} not found in notebook code"
            )

    def test_uns_key_names_in_code(self) -> None:
        """Check that doublet_contract uns key and its sub-keys appear."""
        all_sources = _all_code_sources()
        all_code = "\n".join(all_sources)
        assert 'uns["doublet_contract"]' in all_code, (
            "Required uns key doublet_contract not found"
        )
        for sub_key in ["per_sample_diagnostics", "needs_review", "skip_reason"]:
            assert sub_key in all_code, (
                f"Required uns sub-key {sub_key} not found in notebook code"
            )

    def test_qc_report_fields_in_code(self) -> None:
        qc_report_cell = _cell("# === QC 报告摘要 ===")
        for field in ["n_singlet", "n_uncertain", "n_doublet", "n_excluded"]:
            assert f'"{field}"' in qc_report_cell or f"'{field}'" in qc_report_cell, (
                f"qc_report missing field: {field}"
            )
