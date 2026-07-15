"""P3-b leaf tests for 06_annotated.ipynb: 决策6 (6a-6d) + 全簇PI闸门.

AST/结构断言 verify that the notebook source respects:
- 6a: four provenance annotation fields present
- 6b: LLM suggested never writes final (suggested-only)
- 6c: GREEN confidence alone does NOT auto-produce final; final only in gate cell
- 6d: Cluster_N placeholder only in pi_confirmed/unresolved, never as final
- 红线④: no dir()-based flow control; final_gate_passed initialized in PARAMS
- 红线: no hardcoded old metadata strings in notes cell

行为 exec 断言 construct fake environment and run cells to verify runtime behavior:
- GREEN alone does not produce final
- Full PI confirmation produces final
- Partial confirmation blocks final
- ACCEPT_ALL_SUGGESTED requires explicit PI_CONFIRMED
- Version mismatch blocks final

互斥: 禁 git commit SHA。使用 marker + AST/结构断言。
每个校验都配有"应失败"的真实错误输入（禁恒真/空校验，P2-a 教训）。
"""

import ast
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest  # noqa: F401 — used by pytest fixtures (tmp_path)

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "06_annotated.ipynb"


def _source(cell: dict) -> str:
    return "".join(value) if isinstance(value := cell.get("source", ""), list) else value


def _all_code() -> str:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return "\n".join(
        _source(cell) for cell in notebook["cells"] if cell["cell_type"] == "code"
    )


def _cell(marker: str) -> str:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return next(
        _source(cell)
        for cell in notebook["cells"]
        if cell["cell_type"] == "code" and marker in _source(cell)
    )


# ============================================================
# AST/结构断言（对 cell 源码）
# ============================================================


class TestLLMSuggestedNeverWritesFinal:
    """6b: LLM suggested cell never writes pi_decisions or cell_type_final directly."""

    def test_verdict_cell_never_assigns_pi_decisions(self):
        src = _cell("_verdict_results")
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    target_str = ast.unparse(target) if hasattr(ast, "unparse") else ast.dump(target)
                    if "pi_decisions" in target_str:
                        value_str = ast.unparse(node.value) if hasattr(ast, "unparse") else ast.dump(node.value)
                        assert "_verdict_results" not in value_str, (
                            f"verdict cell assigns pi_decisions from _verdict_results: {value_str[:120]}"
                        )

    def test_verdict_cell_never_assigns_cell_type_final(self):
        src = _cell("_verdict_results")
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    target_str = ast.unparse(target) if hasattr(ast, "unparse") else ast.dump(target)
                    assert "cell_type_final_" not in target_str, (
                        f"_verdict_results cell assigns cell_type_final: {target_str[:120]}"
                    )

    def test_verdict_cell_does_not_write_final_via_green(self):
        assert "cell_type_final" not in _cell("_verdict_results") or "cell_type_final_{" in _cell("_verdict_results"), (
            "cell_type_final should only appear as f-string creation in verdict cell, not assignment"
        )


class TestFinalColumnOnlyInGateCell:
    """6c: cell_type_final only written in the PI-final gate cell."""

    def test_final_assignment_only_in_gate_cell(self):
        notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        final_assign_cells = set()
        for cell in notebook["cells"]:
            if cell["cell_type"] != "code":
                continue
            src = _source(cell)
            if "cell_type_final_" not in src:
                continue
            tree = ast.parse(src)
            for node in ast.walk(tree):
                if isinstance(node, ast.Assign):
                    for target in node.targets:
                        target_str = ast.unparse(target) if hasattr(ast, "unparse") else ast.dump(target)
                        if "cell_type_final_" in target_str:
                            final_assign_cells.add(cell.get("id", "?"))

        assert final_assign_cells.issubset({"gate_p3b_final"}), (
            f"cell_type_final assigned in non-gate cells: {final_assign_cells}"
        )

    def test_final_assignment_under_final_gate_passed_guard(self):
        """The cell_type_final assignment in gate cell is within `if final_gate_passed` block."""
        src = _cell("PI-final 闸门")
        tree = ast.parse(src)

        # Build parent map
        parent_map = {}
        for parent in ast.walk(tree):
            for child in ast.iter_child_nodes(parent):
                parent_map[child] = parent

        # Find all cell_type_final assignments and check they are inside if final_gate_passed
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    target_str = ast.unparse(target) if hasattr(ast, "unparse") else ast.dump(target)
                    if "cell_type_final_" not in target_str:
                        continue
                    # Walk up parents to find an `if final_gate_passed` test
                    current = node
                    found_guard = False
                    while current in parent_map:
                        current = parent_map[current]
                        if isinstance(current, ast.If):
                            test_str = ast.unparse(current.test) if hasattr(ast, "unparse") else ast.dump(current.test)
                            if "final_gate_passed" in test_str:
                                found_guard = True
                                break
                    assert found_guard, (
                        f"cell_type_final assignment not under final_gate_passed guard: "
                        f"{ast.unparse(node) if hasattr(ast, 'unparse') else ast.dump(node)[:200]}"
                    )


class TestFourProvenanceFieldsPresent:
    """6a: four provenance annotation fields present in joined code."""

    def test_five_field_prefixes_present(self):
        joined = _all_code()
        required = [
            "cell_type_marker_",
            "cell_type_llm_suggested_",
            "cell_type_pi_confirmed_",
            "cell_type_final_",
            "annotation_provenance_",
        ]
        for prefix in required:
            assert prefix in joined, f"Missing field prefix in notebook: {prefix}"


class TestClusterNIsPlaceholderNotFinal:
    """6d: Cluster_N placeholder only in pi_confirmed/unresolved context, never as final."""

    def test_cluster_n_never_in_final_assignment(self, tmp_path):
        """6d 红线（行为级）：Cluster_N 占位符绝不能进入 final cell type 列。

        AST 搜索字面量 `cell_type_final_` 在此处不可靠——gate cell 实际用变量
        `_final_col`（`adata.obs[_final_col] = ...`）赋值，AST 展开后目标从不含
        该字面量，导致校验恒过（空跑）。改为行为级 exec：真正运行 decision + gate
        cell，断言未确认簇不会把 Cluster_N 占位写进 final 列。

        场景一（某簇 unresolved）：final 列根本不应被创建。
        场景二（全簇 resolved + PI 确认）：final 列存在，且取值里没有任何
        以 `Cluster_` 前缀的占位标签泄漏。
        """
        # 场景一：2 簇但只确认 1 簇（簇 1 未确认 → decision_source='unresolved'）
        env, adata = _build_gate_env(tmp_path, pi_confirmed=True,
                                     decisions={"0": "T_cell"})
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)
        final_col = f"cell_type_final_{env['ANNOTATION_OUTPUT_VERSION']}"
        assert env["final_gate_passed"] is False, (
            "存在未确认簇时 final_gate_passed 必须为 False"
        )
        assert final_col not in adata.obs.columns, (
            f"未确认簇存在时 {final_col} 不应被创建（Cluster_N 占位不得落 final）"
        )

        # 场景二：全簇确认 + PI_CONFIRMED=True → final 列生成且无 Cluster_ 占位泄漏
        env2, adata2 = _build_gate_env(tmp_path, pi_confirmed=True,
                                       decisions={"0": "T_cell", "1": "B_cell"})
        _exec_cell("=== PI 决策输入", env2)
        _exec_cell("PI-final 闸门", env2)
        assert env2["final_gate_passed"] is True, (
            "全簇确认 + PI_CONFIRMED=True 应通过闸门"
        )
        assert final_col in adata2.obs.columns, (
            f"闸门通过后 {final_col} 应存在: columns={list(adata2.obs.columns)}"
        )
        final_values = [str(v) for v in adata2.obs[final_col].unique()]
        leaked = [v for v in final_values if v.startswith("Cluster_")]
        assert not leaked, (
            f"Cluster_N 占位符泄漏进 final 列: {leaked}（final 取值={final_values}）"
        )

    def test_cluster_n_appears_in_placeholder_context(self):
        joined = _all_code()
        assert "Cluster_" in joined, "Cluster_N placeholder expected but not found"


class TestNoDirBasedFlow:
    """红线④: no dir()-based flow control; final_gate_passed initialized in PARAMS."""

    def test_no_dir_in_condition(self):
        joined = _all_code()
        tree = ast.parse(joined)
        for node in ast.walk(tree):
            if isinstance(node, ast.Compare):
                left_str = ast.unparse(node.left) if hasattr(ast, "unparse") else ast.dump(node.left)
                if "dir()" in left_str:
                    for op in node.ops:
                        if isinstance(op, ast.In):
                            raise AssertionError(
                                f"dir()-based flow control found: "
                                f"{ast.unparse(node) if hasattr(ast, 'unparse') else ast.dump(node)[:200]}"
                            )

    def test_final_gate_passed_initialized_in_params(self):
        params = _cell("# === PARAMS ===")
        assert "final_gate_passed = False" in params, "final_gate_passed not initialized to False in PARAMS"

    def test_pi_confirmed_initialized_in_params(self):
        params = _cell("# === PARAMS ===")
        assert "PI_CONFIRMED" in params, "PI_CONFIRMED not found in PARAMS"


class TestMetadataHasNoHardcodedOldStrings:
    """Notes cell (ee56c5ca) has no hardcoded old metadata strings."""

    def test_no_old_strings_in_notes_cell(self):
        notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        notes_cell = next(c for c in notebook["cells"] if c.get("id") == "ee56c5ca")
        src = _source(notes_cell)
        assert '"PI manual review"' not in src, "Old hardcoded string in notes cell"
        assert "PI manual review of LLM verdicts" not in src, "Old hardcoded string in notes cell"
        assert "final labels reflect domain expertise" not in src, "Old hardcoded string in notes cell"

    def test_notes_cell_references_gate_status(self):
        notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        notes_cell = next(c for c in notebook["cells"] if c.get("id") == "ee56c5ca")
        src = _source(notes_cell)
        has_reference = ("final_gate_passed" in src or "annotation_gate" in src
                         or "confirmation" in src.lower() or "unresolved" in src.lower())
        assert has_reference, "Notes cell should reference gate status or have conditional logic"


class TestVersionFieldsDecoupled:
    """PARAMS has ANNOTATION_OUTPUT_VERSION and UPSTREAM_LABEL_VERSION, no bare OUTPUT_VERSION =."""

    def test_version_fields_are_decoupled(self):
        params = _cell("# === PARAMS ===")
        assert "ANNOTATION_OUTPUT_VERSION" in params, "ANNOTATION_OUTPUT_VERSION missing from PARAMS"
        assert "UPSTREAM_LABEL_VERSION" in params, "UPSTREAM_LABEL_VERSION missing from PARAMS"
        lines = params.split("\n")
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("OUTPUT_VERSION") and "=" in stripped:
                if not stripped.startswith("#"):
                    raise AssertionError(
                        f"Bare OUTPUT_VERSION = assignment found in PARAMS: {stripped[:80]}"
                    )


# ============================================================
# 行为 exec 断言（构造 fake env 跑 gate/decision cell）
# ============================================================


def _safe_sort_key(x):
    try:
        return (0, int(x))
    except (ValueError, TypeError):
        return (1, str(x))


class _FakeAdata:
    """Minimal adata mock for gate cell behavior tests."""

    def __init__(self, n_clusters=2):
        self.X = None
        n_cells = n_clusters * 2
        self.n_obs = n_cells
        self.n_vars = 2000
        clusters = [str(i) for i in range(n_clusters) for _ in range(2)]
        self.obs = pd.DataFrame({"leiden_res_0.6": clusters})
        self.obsm = {}
        self.uns = {}

    @property
    def shape(self):
        return (self.n_obs, self.n_vars)


def _write_pi_csv(csv_path: Path, decisions: dict, version: str = "v1",
                  all_decisions: dict = None):
    """Write a PI confirmation CSV file that the decision cell can read.

    Args:
        csv_path: path to write CSV
        decisions: {cluster: label} for resolved clusters
        version: annotation_version column value
        all_decisions: optional {cluster: (label, decision_type)} for explicit decision types
    """
    clusters = set()
    if all_decisions:
        clusters.update(all_decisions.keys())
    if decisions:
        clusters.update(decisions.keys())
    # Ensure clusters 0..N-1 are all present if decisions cover them
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cluster", "suggested_label", "pi_final_label",
                         "decision", "note", "annotation_version"])
        for cl in sorted(clusters, key=_safe_sort_key):
            if all_decisions and cl in all_decisions:
                label, dec = all_decisions[cl]
            elif cl in decisions:
                label, dec = decisions[cl], "modify"
            else:
                label, dec = "", "unresolved"
            writer.writerow([cl, "", label, dec, "", version])


def _build_gate_env(tmp_path: Path, pi_confirmed=False, decisions=None,
                    csv_version="v1", accept_all=False,
                    verdict_results=None):
    """Build an exec environment for gate + decision cells.

    Decisions are written to a temp CSV file and loaded via PI_CONFIRMATION_CSV.
    """
    n_clusters = max(len(decisions), 2) if decisions else 2
    adata = _FakeAdata(n_clusters=n_clusters)

    # Build _verdict_results
    if verdict_results is None:
        _verdict_results = {}
        for i in range(n_clusters):
            _verdict_results[str(i)] = {
                "label": f"CellType_{i}",
                "confidence": "HIGH",
                "reasoning": f"Test reasoning for cluster {i}",
                "conflict_analysis": "一致",
                "suggestion": "可直接接受",
            }
    else:
        _verdict_results = verdict_results

    # Write temp CSV if decisions provided
    csv_path = ""
    if decisions is not None:
        csv_path = str(tmp_path / "pi_confirm.csv")
        _write_pi_csv(Path(csv_path), decisions, version=csv_version)

    env = {
        "os": __import__("os"),
        "np": np,
        "pd": pd,
        "sp": __import__("scipy.sparse"),
        "Path": Path,
        "json": __import__("json"),
        "datetime": __import__("datetime"),
        "hashlib": __import__("hashlib"),
        "adata": adata,
        "LEIDEN_COL": "leiden_res_0.6",
        "ANNOTATION_OUTPUT_VERSION": "v1",
        "PI_CONFIRMED": pi_confirmed,
        "ACCEPT_ALL_SUGGESTED": accept_all,
        "final_gate_passed": False,
        "PI_CONFIRMATION_CSV": csv_path,
        "_verdict_results": _verdict_results,
        "_pi_csv_version": csv_version,
        "_safe_sort_key": _safe_sort_key,
        "_label_cols": [],
    }
    return env, adata


def _exec_cell(marker: str, env: dict) -> None:
    src = _cell(marker)
    exec(src, env)


class TestGreenAloneDoesNotProduceFinal:
    """6c: PI_CONFIRMED=False + all GREEN -> final_gate_passed is False."""

    def test_green_alone_no_final(self, tmp_path):
        env, adata = _build_gate_env(tmp_path, pi_confirmed=False, decisions={})
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)

        assert not env["final_gate_passed"], (
            f"final_gate_passed should be False when PI_CONFIRMED=False, got {env['final_gate_passed']}"
        )
        final_col = f"cell_type_final_{env['ANNOTATION_OUTPUT_VERSION']}"
        assert final_col not in adata.obs.columns, (
            f"{final_col} should NOT exist when gate is not passed"
        )

    def test_green_with_pi_confirmed_false_rejected(self, tmp_path):
        env, adata = _build_gate_env(tmp_path, pi_confirmed=False,
                                     decisions={"0": "T_cell", "1": "B_cell"})
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)

        assert not env["final_gate_passed"], (
            "final_gate_passed should be False: PI_CONFIRMED=False even with all decisions filled"
        )


class TestFullPIConfirmationProducesFinal:
    """PI_CONFIRMED=True + all clusters covered + no unresolved -> final_gate_passed=True."""

    def test_full_confirmation_produces_final(self, tmp_path):
        env, adata = _build_gate_env(tmp_path, pi_confirmed=True,
                                     decisions={"0": "T_cell", "1": "B_cell"},
                                     csv_version="v1")
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)

        assert env["final_gate_passed"], (
            "final_gate_passed should be True when PI_CONFIRMED=True and all clusters resolved"
        )
        final_col = f"cell_type_final_{env['ANNOTATION_OUTPUT_VERSION']}"
        assert final_col in adata.obs.columns, (
            f"{final_col} should exist when gate is passed: columns={list(adata.obs.columns)}"
        )
        unique_final = set(adata.obs[final_col].unique())
        assert unique_final.issubset({"T_cell", "B_cell"}), (
            f"final column should contain PI decisions, got {unique_final}"
        )

    def test_full_confirmation_notes_record_confirmation(self, tmp_path):
        env, adata = _build_gate_env(tmp_path, pi_confirmed=True,
                                     decisions={"0": "T_cell", "1": "B_cell"})
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)

        notes_key = f"cell_type_final_{env['ANNOTATION_OUTPUT_VERSION']}_notes"
        assert notes_key in adata.uns, f"Notes key {notes_key} not in adata.uns"
        notes = adata.uns[notes_key]
        assert notes.get("status") == "final", f"Notes status should be 'final', got {notes.get('status')}"
        assert notes.get("annotation_version") == env["ANNOTATION_OUTPUT_VERSION"]
        assert notes.get("all_clusters_resolved") is True


class TestPartialConfirmationBlocksFinal:
    """PI_CONFIRMED=True but partial resolution -> final_gate_passed=False."""

    def test_partial_confirmation_blocks_final(self, tmp_path):
        # 3 clusters, only 2 decided
        env, adata = _build_gate_env(tmp_path, pi_confirmed=True,
                                     decisions={"0": "T_cell", "1": "B_cell"})
        # Add a 3rd cluster
        adata.obs = pd.DataFrame({
            "leiden_res_0.6": ["0", "0", "1", "1", "2", "2"]
        })
        adata.n_obs = 6
        env["_verdict_results"]["2"] = {
            "label": "Myeloid", "confidence": "HIGH",
            "reasoning": "Test", "conflict_analysis": "一致", "suggestion": "接受",
        }
        # Rewrite CSV with only 2 decisions out of 3 clusters
        csv_path = str(tmp_path / "pi_confirm.csv")
        _write_pi_csv(Path(csv_path), {"0": "T_cell", "1": "B_cell"})
        env["PI_CONFIRMATION_CSV"] = csv_path

        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)

        assert not env["final_gate_passed"], (
            "final_gate_passed should be False with unresolved cluster 2"
        )
        unresolved = adata.uns.get("annotation_gate", {}).get("unresolved_clusters", [])
        assert "2" in unresolved, f"Cluster 2 should be unresolved, got {unresolved}"
        final_col = f"cell_type_final_{env['ANNOTATION_OUTPUT_VERSION']}"
        assert final_col not in adata.obs.columns, "cell_type_final should not exist"


class TestAcceptAllRequiresExplicitConfirm:
    """ACCEPT_ALL_SUGGESTED only works when PI_CONFIRMED=True."""

    def test_accept_all_without_pi_confirmed_blocked(self, tmp_path):
        env, adata = _build_gate_env(tmp_path, pi_confirmed=False,
                                     decisions={},
                                     accept_all=True)
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)

        assert not env["final_gate_passed"], (
            "ACCEPT_ALL_SUGGESTED=True should not bypass PI_CONFIRMED=False"
        )

    def test_accept_all_with_pi_confirmed_works(self, tmp_path):
        env, adata = _build_gate_env(tmp_path, pi_confirmed=True,
                                     decisions={},
                                     accept_all=True)
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)

        assert env["final_gate_passed"], (
            "ACCEPT_ALL_SUGGESTED=True + PI_CONFIRMED=True should pass the gate"
        )
        prov_key = f"annotation_provenance_{env['ANNOTATION_OUTPUT_VERSION']}"
        prov = adata.uns.get(prov_key, {})
        for cl in ["0", "1"]:
            assert prov.get(cl, {}).get("decision_source") == "accept", (
                f"Cluster {cl} decision_source should be 'accept', got {prov.get(cl, {})}"
            )


class TestVersionMismatchBlocksFinal:
    """Version inconsistency between CSV and notebook blocks final."""

    def test_version_mismatch_blocks_final(self, tmp_path):
        env, adata = _build_gate_env(tmp_path, pi_confirmed=True,
                                     decisions={"0": "T_cell", "1": "B_cell"},
                                     csv_version="v2")  # mismatched: notebook is v1
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)

        assert not env["final_gate_passed"], (
            "Version mismatch (CSV=v2, notebook=v1) should block final"
        )

    def test_version_match_with_csv_passes(self, tmp_path):
        env, adata = _build_gate_env(tmp_path, pi_confirmed=True,
                                     decisions={"0": "T_cell", "1": "B_cell"},
                                     csv_version="v1")  # matched
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)

        assert env["final_gate_passed"], (
            "Version match (CSV=v1, notebook=v1) should pass"
        )


class TestChecksAreNotVacuousTautologies:
    """P2-a 教训：每个校验必须对真实错误输入返回 False（禁恒真/空校验）."""

    def test_green_check_fails_with_wrong_input(self, tmp_path):
        env, adata = _build_gate_env(tmp_path, pi_confirmed=False, decisions={})
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)
        assert env["final_gate_passed"] is False, (
            "PI_CONFIRMED=False must make gate False (NOT vacuous)"
        )

    def test_version_check_fails_with_wrong_input(self, tmp_path):
        env, adata = _build_gate_env(tmp_path, pi_confirmed=True,
                                     decisions={"0": "T_cell", "1": "B_cell"},
                                     csv_version="v99")
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)
        assert env["final_gate_passed"] is False, (
            "v99 vs v1 version mismatch must fail (NOT vacuous)"
        )

    def test_full_gate_pass_actually_passes(self, tmp_path):
        env, adata = _build_gate_env(tmp_path, pi_confirmed=True,
                                     decisions={"0": "T_cell", "1": "B_cell"})
        _exec_cell("=== PI 决策输入", env)
        _exec_cell("PI-final 闸门", env)
        assert env["final_gate_passed"] is True, (
            "Correct inputs (PI_CONFIRMED=True, all resolved, version match) must pass"
        )


class TestCheckpointCellReadsFinalGatePassed:
    """Verify checkpoint cell uses final_gate_passed for needs_review."""

    def test_checkpoint_needs_review_reads_final_gate_passed(self):
        src = _cell("Stage 06 draft checkpoint")
        assert "needs_review = not final_gate_passed" in src, (
            "checkpoint cell must use 'needs_review = not final_gate_passed'"
        )

    def test_checkpoint_has_annotation_gate_in_payload(self):
        src = _cell("Stage 06 draft checkpoint")
        assert '"annotation_gate"' in src, (
            "checkpoint cell must include annotation_gate in manifest_payload"
        )
