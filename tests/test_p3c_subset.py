# ruff: noqa: E401,E501,E701,E702,E731,I001
"""P3-c leaf tests: subset annotation gate + reflow gate (behavioral exec, no AST stub)."""

import ast
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scrna_integration.run_contract import (
    atomic_write_json,
    prepare_run,
    sha256_file,
    validate_artifacts,
    validate_checkpoint,
)

NB = Path(__file__).parents[1] / "notebooks/06c_subset.ipynb"


def _read_nb():
    return json.loads(NB.read_text(encoding="utf-8"))


def _code_cells():
    return [
        "".join(c.get("source", []))
        for c in _read_nb()["cells"]
        if c["cell_type"] == "code"
    ]


def _all_source():
    return "\n\n".join(_code_cells())


def _fn(name, ns):
    tree = ast.parse(_all_source())
    node = next(
        n for n in tree.body
        if isinstance(n, ast.FunctionDef) and n.name == name
    )
    exec(compile(ast.Module([node], []), str(NB), "exec"), ns)
    return ns[name]


def _make_adata(n_cells=5, obs_dict=None, uns_dict=None):
    x = sparse.csr_matrix(np.eye(n_cells, 10, dtype=np.float32))
    obs = pd.DataFrame(obs_dict or {}, index=[f"c{i}" for i in range(n_cells)])
    a = ad.AnnData(x, obs=obs)
    if uns_dict:
        a.uns.update(uns_dict)
    return a


# ---------------------------------------------------------------------------
# 1. subset gate blocked — PI not confirmed
# ---------------------------------------------------------------------------

def test_subset_gate_blocked_pi_not_confirmed():
    """SUBSET_PI_CONFIRMED=False with all pi_decisions → gate returns False, no final col."""
    code_cells = _code_cells()
    gate_idx = None
    for i, src in enumerate(code_cells):
        if "annotation_gate_subset" in src and "subset_final_gate_passed" in src:
            gate_idx = i
            break
    if gate_idx is None:
        pytest.skip("annotation gate cell not found")

    ns = {
        "pd": pd, "np": np, "os": __import__("os"),
        "datetime": __import__("datetime"), "hashlib": __import__("hashlib"),
        "json": json, "Path": Path, "re": re,
        "UPSTREAM_LABEL_VERSION": "v1", "SUBSET_OUTPUT_VERSION": "v1",
        "MAIN_OUTPUT_VERSION": "v2", "UPSTREAM_LABEL_COL": "cell_type_final_v1",
        "SUBSET_LABEL_COL": "cell_type_final_subset_v1",
        "MAIN_LABEL_COL": "cell_type_unified_v2",
        "LEIDEN_COL": "leiden_res_0.6",
        "SUBSET_PI_CONFIRMED": False,
        "MAIN_REFLOW_CONFIRMED": False,
        "SUBSET_ACCEPT_ALL_SUGGESTED": False,
        "subset_final_gate_passed": False,
        "main_reflow_gate_passed": False,
        "SUBSET_PI_CONFIRMATION_CSV": "",
        "SUBSET_MARKER_COL": "cell_type_marker_subset_v1",
        "SUBSET_LLM_COL": "cell_type_llm_suggested_subset_v1",
        "SUBSET_PI_CONFIRMED_COL": "cell_type_pi_confirmed_subset_v1",
        "SUBSET_PROV_KEY": "annotation_provenance_subset_v1",
        "_safe_sort_key": lambda x: (0, int(x)) if str(x).lstrip("-").isdigit() else (1, str(x)),
        "_pi_csv_version": "v1",
    }

    subset = _make_adata(
        6,
        {
            "leiden_res_0.6": ["0", "0", "1", "1", "2", "2"],
            "cell_type_pi_confirmed_subset_v1": [
                "CD4_Tcm", "CD4_Tcm", "CD8_Tem", "CD8_Tem", "Treg", "Treg",
            ],
        },
        {
            "annotation_provenance_subset_v1": {
                "0": {"pi_decision": "CD4_Tcm", "decision_source": "accept"},
                "1": {"pi_decision": "CD8_Tem", "decision_source": "accept"},
                "2": {"pi_decision": "Treg", "decision_source": "accept"},
            }
        },
    )
    subset.obs["leiden_res_0.6"] = subset.obs["leiden_res_0.6"].astype("category")
    ns["adata_sub"] = subset

    exec(compile(code_cells[gate_idx], f"<gate_{gate_idx}>", "exec"), ns)

    assert ns["subset_final_gate_passed"] is False
    assert "cell_type_final_subset_v1" not in subset.obs.columns
    notes_key = "cell_type_final_subset_v1_notes"
    assert subset.uns.get(notes_key, {}).get("status") == "draft-suggested-only"


# ---------------------------------------------------------------------------
# 2. subset gate blocked — unresolved clusters
# ---------------------------------------------------------------------------

def test_subset_gate_blocked_unresolved():
    """PI confirmed but some clusters unresolved → gate blocked."""
    code_cells = _code_cells()
    gate_idx = None
    for i, src in enumerate(code_cells):
        if "annotation_gate_subset" in src and "subset_final_gate_passed" in src:
            gate_idx = i
            break
    if gate_idx is None:
        pytest.skip("annotation gate cell not found")

    ns = {
        "pd": pd, "np": np, "os": __import__("os"),
        "datetime": __import__("datetime"), "hashlib": __import__("hashlib"),
        "json": json, "Path": Path, "re": re,
        "UPSTREAM_LABEL_VERSION": "v1", "SUBSET_OUTPUT_VERSION": "v1",
        "MAIN_OUTPUT_VERSION": "v2", "UPSTREAM_LABEL_COL": "cell_type_final_v1",
        "SUBSET_LABEL_COL": "cell_type_final_subset_v1",
        "MAIN_LABEL_COL": "cell_type_unified_v2",
        "LEIDEN_COL": "leiden_res_0.6",
        "SUBSET_PI_CONFIRMED": True,
        "MAIN_REFLOW_CONFIRMED": False,
        "SUBSET_ACCEPT_ALL_SUGGESTED": False,
        "subset_final_gate_passed": False,
        "main_reflow_gate_passed": False,
        "SUBSET_PI_CONFIRMATION_CSV": "",
        "SUBSET_MARKER_COL": "cell_type_marker_subset_v1",
        "SUBSET_LLM_COL": "cell_type_llm_suggested_subset_v1",
        "SUBSET_PI_CONFIRMED_COL": "cell_type_pi_confirmed_subset_v1",
        "SUBSET_PROV_KEY": "annotation_provenance_subset_v1",
        "_safe_sort_key": lambda x: (0, int(x)) if str(x).lstrip("-").isdigit() else (1, str(x)),
        "_pi_csv_version": "v1",
    }

    # Cluster 2 is unresolved
    subset = _make_adata(
        6,
        {
            "leiden_res_0.6": ["0", "0", "1", "1", "2", "2"],
            "cell_type_pi_confirmed_subset_v1": [
                "CD4_Tcm", "CD4_Tcm", "CD8_Tem", "CD8_Tem", "Cluster_2", "Cluster_2",
            ],
        },
        {
            "annotation_provenance_subset_v1": {
                "0": {"pi_decision": "CD4_Tcm", "decision_source": "accept"},
                "1": {"pi_decision": "CD8_Tem", "decision_source": "accept"},
                "2": {"pi_decision": "", "decision_source": "unresolved"},
            }
        },
    )
    subset.obs["leiden_res_0.6"] = subset.obs["leiden_res_0.6"].astype("category")
    ns["adata_sub"] = subset

    exec(compile(code_cells[gate_idx], f"<gate_{gate_idx}>", "exec"), ns)

    assert ns["subset_final_gate_passed"] is False
    unresolved = subset.uns.get("annotation_gate_subset", {}).get("unresolved_clusters", [])
    assert "2" in unresolved


# ---------------------------------------------------------------------------
# 3. LLM high-confidence does NOT auto-create final (6c semantic)
# ---------------------------------------------------------------------------

def test_llm_suggestions_never_auto_final():
    """LLM suggestions + PI not confirmed → no final column (machine suggestions are advisory)."""
    code_cells = _code_cells()
    pi_idx = None
    for i, src in enumerate(code_cells):
        if "SUBSET_PI_CONFIRMED_COL" in src and "pi_subset_decisions" in src:
            pi_idx = i
            break
    gate_idx = None
    for i, src in enumerate(code_cells):
        if "annotation_gate_subset" in src and "subset_final_gate_passed" in src:
            gate_idx = i
            break
    if pi_idx is None or gate_idx is None:
        pytest.skip("PI decisions or gate cell not found")

    ns = {
        "pd": pd, "np": np, "os": __import__("os"),
        "json": json, "Path": Path, "re": re,
        "UPSTREAM_LABEL_VERSION": "v1", "SUBSET_OUTPUT_VERSION": "v1",
        "MAIN_OUTPUT_VERSION": "v2", "UPSTREAM_LABEL_COL": "cell_type_final_v1",
        "SUBSET_LABEL_COL": "cell_type_final_subset_v1",
        "MAIN_LABEL_COL": "cell_type_unified_v2",
        "LEIDEN_COL": "leiden_res_0.6",
        "SUBSET_PI_CONFIRMED": False,
        "MAIN_REFLOW_CONFIRMED": False,
        "SUBSET_ACCEPT_ALL_SUGGESTED": False,
        "subset_final_gate_passed": False,
        "main_reflow_gate_passed": False,
        "SUBSET_PI_CONFIRMATION_CSV": "",
        "SUBSET_MARKER_COL": "cell_type_marker_subset_v1",
        "SUBSET_LLM_COL": "cell_type_llm_suggested_subset_v1",
        "SUBSET_PI_CONFIRMED_COL": "cell_type_pi_confirmed_subset_v1",
        "SUBSET_PROV_KEY": "annotation_provenance_subset_v1",
        "_safe_sort_key": lambda x: (0, int(x)) if str(x).lstrip("-").isdigit() else (1, str(x)),
        "pi_subset_decisions": {},
        "_pi_from_csv": {},
        "_pi_decision_types": {},
        "_pi_notes": {},
        "_pi_csv_version": "",
        "ANNOTATION_OUTPUT_VERSION": "v1",
        "_verdict_results": {},
        "PI_CONFIRMED": False,
        "ACCEPT_ALL_SUGGESTED": False,
    }

    subset = _make_adata(
        6,
        {
            "leiden_res_0.6": ["0", "0", "1", "1", "2", "2"],
            "cell_type_llm_suggested_subset_v1": [
                "CD4_Tcm", "CD4_Tcm", "CD8_Tem", "CD8_Tem", "Treg", "Treg",
            ],
        },
        {},
    )
    subset.obs["leiden_res_0.6"] = subset.obs["leiden_res_0.6"].astype("category")
    ns["adata_sub"] = subset

    exec(compile(code_cells[pi_idx], f"<pi_{pi_idx}>", "exec"), ns)
    exec(compile(code_cells[gate_idx], f"<gate_{gate_idx}>", "exec"), ns)

    assert ns["subset_final_gate_passed"] is False
    assert "cell_type_final_subset_v1" not in subset.obs.columns


# ---------------------------------------------------------------------------
# 4. Cluster_N placeholder does NOT count as valid label
# ---------------------------------------------------------------------------

def test_cluster_n_placeholder_blocks_gate():
    """pi_confirmed contains Cluster_N → decision_source unresolved → gate blocked."""
    code_cells = _code_cells()
    gate_idx = None
    for i, src in enumerate(code_cells):
        if "annotation_gate_subset" in src and "subset_final_gate_passed" in src:
            gate_idx = i
            break
    if gate_idx is None:
        pytest.skip("annotation gate cell not found")

    ns = {
        "pd": pd, "np": np, "os": __import__("os"),
        "datetime": __import__("datetime"), "hashlib": __import__("hashlib"),
        "json": json, "Path": Path, "re": re,
        "UPSTREAM_LABEL_VERSION": "v1", "SUBSET_OUTPUT_VERSION": "v1",
        "MAIN_OUTPUT_VERSION": "v2", "UPSTREAM_LABEL_COL": "cell_type_final_v1",
        "SUBSET_LABEL_COL": "cell_type_final_subset_v1",
        "MAIN_LABEL_COL": "cell_type_unified_v2",
        "LEIDEN_COL": "leiden_res_0.6",
        "SUBSET_PI_CONFIRMED": True,
        "MAIN_REFLOW_CONFIRMED": False,
        "SUBSET_ACCEPT_ALL_SUGGESTED": False,
        "subset_final_gate_passed": False,
        "main_reflow_gate_passed": False,
        "SUBSET_PI_CONFIRMATION_CSV": "",
        "SUBSET_MARKER_COL": "cell_type_marker_subset_v1",
        "SUBSET_LLM_COL": "cell_type_llm_suggested_subset_v1",
        "SUBSET_PI_CONFIRMED_COL": "cell_type_pi_confirmed_subset_v1",
        "SUBSET_PROV_KEY": "annotation_provenance_subset_v1",
        "_safe_sort_key": lambda x: (0, int(x)) if str(x).lstrip("-").isdigit() else (1, str(x)),
        "_pi_csv_version": "v1",
    }

    subset = _make_adata(
        4,
        {
            "leiden_res_0.6": ["0", "0", "1", "1"],
            "cell_type_pi_confirmed_subset_v1": [
                "CD4_Tcm", "CD4_Tcm", "Cluster_1", "Cluster_1",
            ],
        },
        {
            "annotation_provenance_subset_v1": {
                "0": {"pi_decision": "CD4_Tcm", "decision_source": "accept"},
                "1": {"pi_decision": "", "decision_source": "unresolved"},
            }
        },
    )
    subset.obs["leiden_res_0.6"] = subset.obs["leiden_res_0.6"].astype("category")
    ns["adata_sub"] = subset

    exec(compile(code_cells[gate_idx], f"<gate_{gate_idx}>", "exec"), ns)

    assert ns["subset_final_gate_passed"] is False
    unresolved = subset.uns.get("annotation_gate_subset", {}).get("unresolved_clusters", [])
    assert "1" in unresolved


# ---------------------------------------------------------------------------
# 5. reflow gate blocked — subset gate not passed
# ---------------------------------------------------------------------------

def test_reflow_gate_blocked_subset_not_passed():
    """subset_final_gate_passed=False → main_reflow_gate_passed False."""
    code_cells = _code_cells()
    reflow_idx = None
    for i, src in enumerate(code_cells):
        if ("main_reflow_gate_passed" in src and "MAIN_REFLOW_CONFIRMED" in src
                and "subset_final_gate_passed" in src and "MAIN_OUTPUT_VERSION" in src):
            # Check it's the reflow gate, not the output cell
            if "= bool(" in src or "main_reflow_gate_passed =" in src:
                reflow_idx = i
                break
    if reflow_idx is None:
        pytest.skip("reflow gate cell not found")

    ns = {
        "subset_final_gate_passed": False,
        "MAIN_REFLOW_CONFIRMED": True,
        "MAIN_OUTPUT_VERSION": "v2",
        "UPSTREAM_LABEL_VERSION": "v1",
        "main_reflow_gate_passed": False,
    }

    exec(compile(code_cells[reflow_idx], f"<reflow_{reflow_idx}>", "exec"), ns)

    assert ns["main_reflow_gate_passed"] is False


# ---------------------------------------------------------------------------
# 6. reflow gate blocked — MAIN_REFLOW_CONFIRMED=False
# ---------------------------------------------------------------------------

def test_reflow_gate_blocked_pi_not_confirmed():
    """subset gate passed but MAIN_REFLOW_CONFIRMED=False → main_reflow_gate_passed False."""
    code_cells = _code_cells()
    reflow_idx = None
    for i, src in enumerate(code_cells):
        if ("main_reflow_gate_passed" in src and "MAIN_REFLOW_CONFIRMED" in src
                and "subset_final_gate_passed" in src and "MAIN_OUTPUT_VERSION" in src):
            if "= bool(" in src or "main_reflow_gate_passed =" in src:
                reflow_idx = i
                break
    if reflow_idx is None:
        pytest.skip("reflow gate cell not found")

    ns = {
        "subset_final_gate_passed": True,
        "MAIN_REFLOW_CONFIRMED": False,
        "MAIN_OUTPUT_VERSION": "v2",
        "UPSTREAM_LABEL_VERSION": "v1",
        "main_reflow_gate_passed": False,
    }

    exec(compile(code_cells[reflow_idx], f"<reflow_{reflow_idx}>", "exec"), ns)

    assert ns["main_reflow_gate_passed"] is False


# ---------------------------------------------------------------------------
# 7. version decoupling — MAIN_OUTPUT_VERSION==UPSTREAM_LABEL_VERSION blocks reflow
# ---------------------------------------------------------------------------

def test_reflow_gate_blocked_version_not_bumped():
    """MAIN_OUTPUT_VERSION==UPSTREAM_LABEL_VERSION → main_reflow_gate_passed False."""
    code_cells = _code_cells()
    reflow_idx = None
    for i, src in enumerate(code_cells):
        if ("main_reflow_gate_passed" in src and "MAIN_REFLOW_CONFIRMED" in src
                and "subset_final_gate_passed" in src and "MAIN_OUTPUT_VERSION" in src):
            if "= bool(" in src or "main_reflow_gate_passed =" in src:
                reflow_idx = i
                break
    if reflow_idx is None:
        pytest.skip("reflow gate cell not found")

    ns = {
        "subset_final_gate_passed": True,
        "MAIN_REFLOW_CONFIRMED": True,
        "MAIN_OUTPUT_VERSION": "v1",
        "UPSTREAM_LABEL_VERSION": "v1",
        "main_reflow_gate_passed": False,
    }

    exec(compile(code_cells[reflow_idx], f"<reflow_{reflow_idx}>", "exec"), ns)

    assert ns["main_reflow_gate_passed"] is False


# ---------------------------------------------------------------------------
# 8. annotation confirmed but reflow not — subset has final, main untouched
# ---------------------------------------------------------------------------

def test_annotation_confirmed_main_not_reflowed():
    """SUBSET_PI_CONFIRMED=True + all resolved → subset final exists. Main not touched."""
    code_cells = _code_cells()
    gate_idx = None
    for i, src in enumerate(code_cells):
        if "annotation_gate_subset" in src and "subset_final_gate_passed" in src:
            gate_idx = i
            break
    if gate_idx is None:
        pytest.skip("annotation gate cell not found")

    ns = {
        "pd": pd, "np": np, "os": __import__("os"),
        "datetime": __import__("datetime"), "hashlib": __import__("hashlib"),
        "json": json, "Path": Path, "re": re,
        "UPSTREAM_LABEL_VERSION": "v1", "SUBSET_OUTPUT_VERSION": "v1",
        "MAIN_OUTPUT_VERSION": "v2", "UPSTREAM_LABEL_COL": "cell_type_final_v1",
        "SUBSET_LABEL_COL": "cell_type_final_subset_v1",
        "MAIN_LABEL_COL": "cell_type_unified_v2",
        "LEIDEN_COL": "leiden_res_0.6",
        "SUBSET_PI_CONFIRMED": True,
        "MAIN_REFLOW_CONFIRMED": False,
        "SUBSET_ACCEPT_ALL_SUGGESTED": False,
        "subset_final_gate_passed": False,
        "main_reflow_gate_passed": False,
        "SUBSET_PI_CONFIRMATION_CSV": "",
        "SUBSET_MARKER_COL": "cell_type_marker_subset_v1",
        "SUBSET_LLM_COL": "cell_type_llm_suggested_subset_v1",
        "SUBSET_PI_CONFIRMED_COL": "cell_type_pi_confirmed_subset_v1",
        "SUBSET_PROV_KEY": "annotation_provenance_subset_v1",
        "_safe_sort_key": lambda x: (0, int(x)) if str(x).lstrip("-").isdigit() else (1, str(x)),
        "_pi_csv_version": "v1",
    }

    subset = _make_adata(
        4,
        {
            "leiden_res_0.6": ["0", "0", "1", "1"],
            "cell_type_pi_confirmed_subset_v1": [
                "CD4_Tcm", "CD4_Tcm", "CD8_Tem", "CD8_Tem",
            ],
        },
        {
            "annotation_provenance_subset_v1": {
                "0": {"pi_decision": "CD4_Tcm", "decision_source": "accept"},
                "1": {"pi_decision": "CD8_Tem", "decision_source": "accept"},
            }
        },
    )
    subset.obs["leiden_res_0.6"] = subset.obs["leiden_res_0.6"].astype("category")
    ns["adata_sub"] = subset

    exec(compile(code_cells[gate_idx], f"<gate_{gate_idx}>", "exec"), ns)

    assert ns["subset_final_gate_passed"] is True
    assert "cell_type_final_subset_v1" in subset.obs.columns


# ---------------------------------------------------------------------------
# 9. CSV version mismatch blocks gate
# ---------------------------------------------------------------------------

def test_csv_version_mismatch_blocks_gate():
    """CSV annotation_version != SUBSET_OUTPUT_VERSION → version_ok False → blocked."""
    code_cells = _code_cells()
    gate_idx = None
    for i, src in enumerate(code_cells):
        if "annotation_gate_subset" in src and "subset_final_gate_passed" in src:
            gate_idx = i
            break
    if gate_idx is None:
        pytest.skip("annotation gate cell not found")

    ns = {
        "pd": pd, "np": np, "os": __import__("os"),
        "datetime": __import__("datetime"), "hashlib": __import__("hashlib"),
        "json": json, "Path": Path, "re": re,
        "UPSTREAM_LABEL_VERSION": "v1", "SUBSET_OUTPUT_VERSION": "v2",
        "MAIN_OUTPUT_VERSION": "v3", "UPSTREAM_LABEL_COL": "cell_type_final_v1",
        "SUBSET_LABEL_COL": "cell_type_final_subset_v2",
        "MAIN_LABEL_COL": "cell_type_unified_v3",
        "LEIDEN_COL": "leiden_res_0.6",
        "SUBSET_PI_CONFIRMED": True,
        "MAIN_REFLOW_CONFIRMED": False,
        "SUBSET_ACCEPT_ALL_SUGGESTED": False,
        "subset_final_gate_passed": False,
        "main_reflow_gate_passed": False,
        "SUBSET_PI_CONFIRMATION_CSV": "dummy.csv",
        "SUBSET_MARKER_COL": "cell_type_marker_subset_v2",
        "SUBSET_LLM_COL": "cell_type_llm_suggested_subset_v2",
        "SUBSET_PI_CONFIRMED_COL": "cell_type_pi_confirmed_subset_v2",
        "SUBSET_PROV_KEY": "annotation_provenance_subset_v2",
        "_safe_sort_key": lambda x: (0, int(x)) if str(x).lstrip("-").isdigit() else (1, str(x)),
        "_pi_csv_version": "v1",
    }

    subset = _make_adata(
        4,
        {
            "leiden_res_0.6": ["0", "0", "1", "1"],
            "cell_type_pi_confirmed_subset_v2": [
                "CD4_Tcm", "CD4_Tcm", "CD8_Tem", "CD8_Tem",
            ],
        },
        {
            "annotation_provenance_subset_v2": {
                "0": {"pi_decision": "CD4_Tcm", "decision_source": "accept"},
                "1": {"pi_decision": "CD8_Tem", "decision_source": "accept"},
            }
        },
    )
    subset.obs["leiden_res_0.6"] = subset.obs["leiden_res_0.6"].astype("category")
    ns["adata_sub"] = subset

    exec(compile(code_cells[gate_idx], f"<gate_{gate_idx}>", "exec"), ns)

    assert ns["subset_final_gate_passed"] is False


# ---------------------------------------------------------------------------
# 10. happy path — both gates pass, full dual-artifact output
# ---------------------------------------------------------------------------

def test_happy_path_both_gates_pass(tmp_path):
    """Both gates pass → cell_type_final_subset created, dual artifacts written."""
    code_cells = _code_cells()
    # Find _write_06c_outputs index
    output_idx = None
    for i, src in enumerate(code_cells):
        if "_write_06c_outputs" in src and "def _write_06c_outputs" in src:
            output_idx = i
            break
    if output_idx is None:
        pytest.skip("_write_06c_outputs function not found")

    ns_base = {
        "Path": Path, "re": re, "pd": pd, "np": np,
        "sc": __import__("scanpy"), "sp": __import__("scipy.sparse"),
        "json": json, "os": __import__("os"), "datetime": __import__("datetime"),
        "prepare_run": prepare_run, "atomic_write_json": atomic_write_json,
        "sha256_file": sha256_file,
        "snapshot_effective_parameters": lambda *a, **k: {},
        "collect_runtime_provenance": lambda *a, **k: {},
        "validate_checkpoint": validate_checkpoint,
        "validate_artifacts": validate_artifacts,
        "_OUTPUT_VERSION_RE": re.compile(r"v[1-9][0-9]*\Z"),
        "_H5AD_BASENAME_RE": re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\.h5ad\Z"),
    }

    ns = dict(ns_base)
    ns.update({
        "RUN_ROOT": tmp_path / "runs", "RUN_ID": "06c-happy",
        "UPSTREAM_RUN_ID": "stage06",
        "UPSTREAM_LABEL_VERSION": "v1", "SUBSET_OUTPUT_VERSION": "v1",
        "MAIN_OUTPUT_VERSION": "v2", "SUBSET_LABELS": ["A"],
        "UPSTREAM_LABEL_COL": "cell_type_final_v1",
        "SUBSET_LABEL_COL": "cell_type_final_subset_v1",
        "MAIN_LABEL_COL": "cell_type_unified_v2",
        "_root": tmp_path,
        "_upstream_manifest_path": tmp_path / "upstream.json",
        "_upstream_manifest_sha256": "0" * 64,
        "_upstream_checkpoint": tmp_path / "upstream.h5ad",
        "_upstream_checkpoint_sha256": "0" * 64,
        "_verified_upstream_index": pd.Index(["c1", "c2", "c3", "c4"]),
        "_verified_upstream_coarse": pd.Series(
            ["A", "A", "B", "B"], index=["c1", "c2", "c3", "c4"]
        ),
        "_verified_upstream_obs_columns": ("cell_type_final_v1",),
    })
    (tmp_path / "upstream.json").write_text("{}")
    (tmp_path / "upstream.h5ad").write_bytes(b"upstream")

    for name in ("_series_equal", "_output_contract", "_manifest_base",
                 "_record_output_failure", "_complete_output_write",
                 "_write_06c_outputs"):
        _fn(name, ns)

    obs = pd.DataFrame(
        {"cell_type_final_v1": ["A", "A", "B", "B"]},
        index=["c1", "c2", "c3", "c4"],
    )
    main = ad.AnnData(sparse.csr_matrix(np.eye(4, dtype=np.float32)), obs=obs.copy())
    main.uns.update({"stage": "06_annotated", "status": "SUCCESS", "run_id": "stage06"})
    subset = main[["c1", "c2"]].copy()
    subset.obs["cell_type_final_subset_v1"] = ["A1", "A2"]
    main.obs["cell_type_final_subset_v1"] = ["A1", "A2", np.nan, np.nan]
    main.obs["cell_type_unified_v2"] = ["A1", "A2", "B", "B"]

    paths, subset_path, main_path = ns["_write_06c_outputs"](subset, main)

    subset_rt = ad.read_h5ad(subset_path)
    main_rt = ad.read_h5ad(main_path)
    assert validate_checkpoint(paths.manifest_path) == subset_path
    artifacts = validate_artifacts(paths.manifest_path)
    assert set(artifacts.keys()) == {"subset", "main_reflow"}
    assert main_rt.obs["cell_type_final_v1"].tolist() == ["A", "A", "B", "B"]
    assert main_rt.obs["cell_type_final_subset_v1"].iloc[:2].tolist() == ["A1", "A2"]
    assert main_rt.obs["cell_type_final_subset_v1"].iloc[2:].isna().all()
    assert main_rt.obs["cell_type_unified_v2"].tolist() == ["A1", "A2", "B", "B"]
    assert set(main_rt.obs.columns) == {
        "cell_type_final_v1", "cell_type_final_subset_v1", "cell_type_unified_v2",
    }
    assert subset_rt.uns["stage"] == "06c_subset"


# ---------------------------------------------------------------------------
# 11. second write failure → atomic rollback retained
# ---------------------------------------------------------------------------

def test_atomic_rollback_second_write_failure(tmp_path):
    """Second write fails → both h5ads unlinked, FAILED manifest, no checkpoint/artifacts."""
    code_cells = _code_cells()
    output_idx = None
    for i, src in enumerate(code_cells):
        if "_write_06c_outputs" in src and "def _write_06c_outputs" in src:
            output_idx = i
            break
    if output_idx is None:
        pytest.skip("_write_06c_outputs not found")

    ns_base = {
        "Path": Path, "re": re, "pd": pd, "np": np,
        "sc": __import__("scanpy"), "sp": __import__("scipy.sparse"),
        "json": json, "os": __import__("os"), "datetime": __import__("datetime"),
        "prepare_run": prepare_run, "atomic_write_json": atomic_write_json,
        "sha256_file": sha256_file,
        "snapshot_effective_parameters": lambda *a, **k: {},
        "collect_runtime_provenance": lambda *a, **k: {},
        "validate_checkpoint": validate_checkpoint,
        "validate_artifacts": validate_artifacts,
        "_OUTPUT_VERSION_RE": re.compile(r"v[1-9][0-9]*\Z"),
        "_H5AD_BASENAME_RE": re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\.h5ad\Z"),
    }
    ns = dict(ns_base)
    ns.update({
        "RUN_ROOT": tmp_path / "runs", "RUN_ID": "06c-rollback",
        "UPSTREAM_RUN_ID": "stage06",
        "UPSTREAM_LABEL_VERSION": "v1", "SUBSET_OUTPUT_VERSION": "v1",
        "MAIN_OUTPUT_VERSION": "v2", "SUBSET_LABELS": ["A"],
        "UPSTREAM_LABEL_COL": "cell_type_final_v1",
        "SUBSET_LABEL_COL": "cell_type_final_subset_v1",
        "MAIN_LABEL_COL": "cell_type_unified_v2",
        "_root": tmp_path,
        "_upstream_manifest_path": tmp_path / "upstream.json",
        "_upstream_manifest_sha256": "0" * 64,
        "_upstream_checkpoint": tmp_path / "upstream.h5ad",
        "_upstream_checkpoint_sha256": "0" * 64,
        "_verified_upstream_index": pd.Index(["c1", "c2", "c3", "c4"]),
        "_verified_upstream_coarse": pd.Series(
            ["A", "A", "B", "B"], index=["c1", "c2", "c3", "c4"]
        ),
        "_verified_upstream_obs_columns": ("cell_type_final_v1",),
    })
    (tmp_path / "upstream.json").write_text("{}")
    (tmp_path / "upstream.h5ad").write_bytes(b"upstream")

    for name in ("_series_equal", "_output_contract", "_manifest_base",
                 "_record_output_failure", "_complete_output_write",
                 "_write_06c_outputs"):
        _fn(name, ns)

    obs = pd.DataFrame(
        {"cell_type_final_v1": ["A", "A", "B", "B"]},
        index=["c1", "c2", "c3", "c4"],
    )
    main = ad.AnnData(sparse.csr_matrix(np.eye(4, dtype=np.float32)), obs=obs.copy())
    main.uns.update({"stage": "06_annotated", "status": "SUCCESS", "run_id": "stage06"})
    subset = main[["c1", "c2"]].copy()
    subset.obs["cell_type_final_subset_v1"] = ["A1", "A2"]
    main.obs["cell_type_final_subset_v1"] = ["A1", "A2", np.nan, np.nan]
    main.obs["cell_type_unified_v2"] = ["A1", "A2", "B", "B"]

    orig_write = main.write_h5ad

    def fail_main(path, **_kw):
        orig_write(path, compression="lzf")
        raise OSError("second write failed")

    main.write_h5ad = fail_main

    with pytest.raises(OSError, match="second write failed"):
        ns["_write_06c_outputs"](subset, main)

    draft = tmp_path / "runs" / "06c-rollback" / "draft"
    m = json.loads((draft / "manifest.json").read_text())
    assert not list(draft.glob("*.h5ad"))
    assert m["stage_status"] == "FAILED"
    assert "checkpoint" not in m
    assert "artifacts" not in m


# ---------------------------------------------------------------------------
# 12. anti-dir-flow-control + PARAMS initial values + no stale literals
# ---------------------------------------------------------------------------

def test_params_have_initial_values():
    """Gate variables defined in PARAMS with explicit initial values, not dir()."""
    code_cells = _code_cells()
    params_cell = code_cells[0]
    for var in ("SUBSET_PI_CONFIRMED", "MAIN_REFLOW_CONFIRMED",
                "subset_final_gate_passed", "main_reflow_gate_passed",
                "SUBSET_ACCEPT_ALL_SUGGESTED", "SUBSET_PI_CONFIRMATION_CSV"):
        assert var in params_cell, f"{var} missing from PARAMS"


def test_no_dir_flow_control_in_gate_branches():
    """Gate cell business branches must not use 'x in dir()' for gate variables."""
    code_cells = _code_cells()
    for src in code_cells:
        if ("subset_final_gate_passed" in src and "annotation_gate_subset" in src) or \
           ("main_reflow_gate_passed" in src and "MAIN_REFLOW_CONFIRMED" in src):
            for line in src.split("\n"):
                s = line.strip()
                if "dir()" in s and not s.startswith("#"):
                    if any(kw in s for kw in ("if ", "elif ")):
                        if any(v in s for v in ("subset_final_gate_passed",
                                                 "main_reflow_gate_passed",
                                                 "SUBSET_PI_CONFIRMED",
                                                 "MAIN_REFLOW_CONFIRMED")):
                            pytest.fail(f"dir() flow control: {s}")


def test_no_hardcoded_column_literals():
    """After setup cell, column names use f-strings, not hardcoded literals."""
    code_cells = _code_cells()
    for src in code_cells[2:]:
        for literal in ("cell_type_final_subset_v1", "cell_type_marker_subset_v1",
                        "cell_type_llm_suggested_subset_v1",
                        "cell_type_pi_confirmed_subset_v1"):
            if literal in src:
                for line in src.split("\n"):
                    if literal in line:
                        s = line.strip()
                        if s.startswith("#"):
                            continue
                        if s.startswith('"') or s.startswith("'"):
                            continue
                        if "=" in s or "obs[" in s:
                            pytest.fail(f"Hardcoded '{literal}': {s[:80]}")


# ---------------------------------------------------------------------------
# 13. _write_06c_subset_draft — draft-only output
# ---------------------------------------------------------------------------

def test_subset_draft_output(tmp_path):
    """_write_06c_subset_draft writes subset only, NEEDS_REVIEW, single artifact."""
    code_cells = _code_cells()
    draft_idx = None
    for i, src in enumerate(code_cells):
        if "_write_06c_subset_draft" in src and "def _write_06c_subset_draft" in src:
            draft_idx = i
            break
    if draft_idx is None:
        pytest.skip("_write_06c_subset_draft not found")

    ns_base = {
        "Path": Path, "re": re, "pd": pd, "np": np,
        "sc": __import__("scanpy"), "sp": __import__("scipy.sparse"),
        "json": json, "os": __import__("os"), "datetime": __import__("datetime"),
        "prepare_run": prepare_run, "atomic_write_json": atomic_write_json,
        "sha256_file": sha256_file,
        "snapshot_effective_parameters": lambda *a, **k: {},
        "collect_runtime_provenance": lambda *a, **k: {},
        "validate_checkpoint": validate_checkpoint,
        "validate_artifacts": validate_artifacts,
        "_OUTPUT_VERSION_RE": re.compile(r"v[1-9][0-9]*\Z"),
        "_H5AD_BASENAME_RE": re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\.h5ad\Z"),
    }
    ns = dict(ns_base)
    ns.update({
        "RUN_ROOT": tmp_path / "runs", "RUN_ID": "06c-draft",
        "UPSTREAM_RUN_ID": "stage06",
        "UPSTREAM_LABEL_VERSION": "v1", "SUBSET_OUTPUT_VERSION": "v1",
        "MAIN_OUTPUT_VERSION": "v2", "SUBSET_LABELS": ["A"],
        "UPSTREAM_LABEL_COL": "cell_type_final_v1",
        "SUBSET_LABEL_COL": "cell_type_final_subset_v1",
        "MAIN_LABEL_COL": "cell_type_unified_v2",
        "SUBSET_PI_CONFIRMED_COL": "cell_type_pi_confirmed_subset_v1",
        "subset_final_gate_passed": False, "main_reflow_gate_passed": False,
        "SUBSET_PI_CONFIRMED": False, "MAIN_REFLOW_CONFIRMED": False,
        "_root": tmp_path,
        "_upstream_manifest_path": tmp_path / "upstream.json",
        "_upstream_manifest_sha256": "0" * 64,
        "_upstream_checkpoint": tmp_path / "upstream.h5ad",
        "_upstream_checkpoint_sha256": "0" * 64,
    })
    (tmp_path / "upstream.json").write_text("{}")
    (tmp_path / "upstream.h5ad").write_bytes(b"upstream")

    _fn("_manifest_base", ns)
    _fn("_write_06c_subset_draft", ns)

    subset = ad.AnnData(
        sparse.csr_matrix(np.eye(3, 10, dtype=np.float32)),
        obs=pd.DataFrame({
            "cell_type_final_subset_v1": ["A1", "A2", np.nan],
            "cell_type_pi_confirmed_subset_v1": ["A1", "A2", "Cluster_2"],
        }, index=["c1", "c2", "c3"]),
    )
    subset.uns["annotation_gate_subset"] = {
        "subset_final_gate_passed": False,
        "unresolved_clusters": ["2"],
    }

    paths, draft_path = ns["_write_06c_subset_draft"](subset)

    assert draft_path.exists()
    manifest = json.loads(paths.manifest_path.read_text())
    assert manifest["stage_status"] == "NEEDS_REVIEW"
    assert len(manifest.get("artifacts", [])) == 1
    assert manifest["artifacts"][0]["role"] == "subset_draft"
    assert manifest["hard_postconditions"]["subset_final_gate_passed"] is False
    assert manifest["gate_status"]["subset_final_gate_passed"] is False
    rt = ad.read_h5ad(draft_path)
    assert rt.uns["stage"] == "06c_subset"
    assert rt.uns["status"] == "NEEDS_REVIEW"


# ---------------------------------------------------------------------------
# 14. _write_06c_subset_draft atomic rollback
# ---------------------------------------------------------------------------

def test_subset_draft_atomic_rollback(tmp_path):
    """_write_06c_subset_draft write failure → unlink + FAILED manifest."""
    code_cells = _code_cells()
    draft_idx = None
    for i, src in enumerate(code_cells):
        if "_write_06c_subset_draft" in src and "def _write_06c_subset_draft" in src:
            draft_idx = i
            break
    if draft_idx is None:
        pytest.skip("_write_06c_subset_draft not found")

    ns_base = {
        "Path": Path, "re": re, "pd": pd, "np": np,
        "sc": __import__("scanpy"), "sp": __import__("scipy.sparse"),
        "json": json, "os": __import__("os"), "datetime": __import__("datetime"),
        "prepare_run": prepare_run, "atomic_write_json": atomic_write_json,
        "sha256_file": sha256_file,
        "snapshot_effective_parameters": lambda *a, **k: {},
        "collect_runtime_provenance": lambda *a, **k: {},
        "validate_checkpoint": validate_checkpoint,
        "validate_artifacts": validate_artifacts,
        "_OUTPUT_VERSION_RE": re.compile(r"v[1-9][0-9]*\Z"),
        "_H5AD_BASENAME_RE": re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\.h5ad\Z"),
    }
    ns = dict(ns_base)
    ns.update({
        "RUN_ROOT": tmp_path / "runs", "RUN_ID": "06c-draft-fail",
        "UPSTREAM_RUN_ID": "stage06",
        "UPSTREAM_LABEL_VERSION": "v1", "SUBSET_OUTPUT_VERSION": "v1",
        "MAIN_OUTPUT_VERSION": "v2", "SUBSET_LABELS": ["A"],
        "UPSTREAM_LABEL_COL": "cell_type_final_v1",
        "SUBSET_LABEL_COL": "cell_type_final_subset_v1",
        "MAIN_LABEL_COL": "cell_type_unified_v2",
        "SUBSET_PI_CONFIRMED_COL": "cell_type_pi_confirmed_subset_v1",
        "subset_final_gate_passed": False, "main_reflow_gate_passed": False,
        "SUBSET_PI_CONFIRMED": False, "MAIN_REFLOW_CONFIRMED": False,
        "_root": tmp_path,
        "_upstream_manifest_path": tmp_path / "upstream.json",
        "_upstream_manifest_sha256": "0" * 64,
        "_upstream_checkpoint": tmp_path / "upstream.h5ad",
        "_upstream_checkpoint_sha256": "0" * 64,
    })
    (tmp_path / "upstream.json").write_text("{}")
    (tmp_path / "upstream.h5ad").write_bytes(b"upstream")

    _fn("_manifest_base", ns)
    _fn("_record_output_failure", ns)
    _fn("_write_06c_subset_draft", ns)

    subset = ad.AnnData(
        sparse.csr_matrix(np.eye(3, 10, dtype=np.float32)),
        obs=pd.DataFrame({
            "cell_type_final_subset_v1": ["A1", "A2", np.nan],
            "cell_type_pi_confirmed_subset_v1": ["A1", "A2", "Cluster_2"],
        }, index=["c1", "c2", "c3"]),
    )
    subset.uns["annotation_gate_subset"] = {"subset_final_gate_passed": False}
    subset.write_h5ad = lambda path, **_: (_ for _ in ()).throw(OSError("write failed"))

    with pytest.raises(OSError, match="write failed"):
        ns["_write_06c_subset_draft"](subset)

    draft = tmp_path / "runs" / "06c-draft-fail" / "draft"
    m = json.loads((draft / "manifest.json").read_text())
    assert m["stage_status"] == "FAILED"
    assert not list(draft.glob("*.h5ad"))
