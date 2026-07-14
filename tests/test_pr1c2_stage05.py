import ast
import json
import re
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import scrna_integration.run_contract as rc

ROOT = Path(__file__).parents[1]
NOTEBOOK = ROOT / "notebooks/05_clustered.ipynb"


def _source(cell: dict) -> str:
    value = cell.get("source", "")
    return "".join(value) if isinstance(value, list) else value


def _cell(marker: str) -> str:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return next(
        _source(cell) for cell in notebook["cells"]
        if cell["cell_type"] == "code" and marker in _source(cell)
    )


class _Adata:
    def __init__(self) -> None:
        self.X = sp.csr_matrix(np.ones((4, 2), dtype=np.float32))
        self.n_obs, self.n_vars = self.X.shape
        self.obs = pd.DataFrame({
            "leiden_res_0.2": ["0", "0", "1", "1"],
            "leiden_res_0.4": ["0", "1", "2", "2"],
        })
        self.obsm = {
            "X_pca": np.ones((4, 2), dtype=np.float32),
            "X_pca_harmony": np.ones((4, 2), dtype=np.float32),
        }
        self.uns: dict = {}
        self.fail_write = False
        self.write_calls = 0

    def write_h5ad(self, path: Path, **_: object) -> None:
        self.write_calls += 1
        Path(path).write_bytes(b"partial" if self.fail_write else b"stage05")
        if self.fail_write:
            raise OSError("disk full")


def _upstream(root: Path) -> tuple[Path, Path]:
    paths = rc.prepare_run(root, "stage04")
    checkpoint = paths.draft_dir / "04.h5ad"
    checkpoint.write_bytes(b"stage04")
    rc.atomic_write_json(paths.manifest_path, {
        "run_id": "stage04", "stage": "04_embedded", "stage_status": "SUCCESS",
        "checkpoint": {"path": checkpoint.name, "sha256": rc.sha256_file(checkpoint)},
    })
    rc.promote_run(paths)
    return paths.promoted_dir / checkpoint.name, paths.promoted_dir / "manifest.json"


def _load_env(tmp_path: Path, case: str = "ok") -> tuple[dict, _Adata]:
    checkpoint, manifest_path = _upstream(tmp_path / "upstream")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if case == "run_id":
        manifest["run_id"] = "other"
    elif case == "stage":
        manifest["stage"] = "03_normalized"
    elif case == "path":
        manifest["checkpoint"]["path"] = "../04.h5ad"
    elif case == "missing_status":
        manifest.pop("stage_status")
    elif case in {"failed", "needs_review"}:
        manifest["stage_status"] = case.upper()
    elif case in {"warnings", "unaccepted_warnings"}:
        manifest["stage_status"] = "SUCCESS_WITH_WARNINGS"
        if case == "warnings":
            manifest["warning_acceptance"] = {
                "accepted_by": "researcher", "accepted_at": "2026-07-13T12:00:00Z"
            }
    rc.atomic_write_json(manifest_path, manifest, overwrite=True)
    if case == "hash":
        checkpoint.write_bytes(b"corrupt")
    adata = _Adata()
    env = {
        "Path": Path,
        "json": json,
        "sp": sp,
        "resume_run": rc.resume_run,
        "sha256_file": rc.sha256_file,
        "validate_checkpoint": rc.validate_checkpoint,
        "sc": SimpleNamespace(read_h5ad=lambda path: adata if Path(path) == checkpoint else None),
    }
    exec(_cell("# === PARAMS ==="), env)
    env.update(UPSTREAM_RUN_ROOT=str(tmp_path / "upstream"), UPSTREAM_RUN_ID="stage04")
    exec(_cell("# === 只读取 verified promoted Stage 04 上游 ==="), env)
    return env, adata


def _final_env(tmp_path: Path) -> tuple[dict, _Adata]:
    env, adata = _load_env(tmp_path)
    env.update(
        Path=Path, re=re, RESOLUTIONS=[0.2, 0.4], RUN_ROOT=str(tmp_path / "runs"),
        RUN_ID="stage05", OUTPUT_FILENAME="05.h5ad", _root=str(tmp_path),
        atomic_write_json=rc.atomic_write_json, collect_runtime_provenance=lambda *_: {"python": "test"},
        determine_stage_status=rc.determine_stage_status, prepare_run=rc.prepare_run,
        snapshot_effective_parameters=rc.snapshot_effective_parameters,
    )
    return env, adata


def test_structure_preserves_sweep_and_defers_selection_and_promotion() -> None:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    code = [_source(cell) for cell in notebook["cells"] if cell["cell_type"] == "code"]
    joined = "\n".join(code)
    ast.parse(joined, filename=str(NOTEBOOK))
    params = _cell("# === PARAMS ===")
    final = _cell("# === Stage 05 draft checkpoint")
    assert all(name in params for name in (
        "UPSTREAM_RUN_ROOT", "UPSTREAM_RUN_ID", "RUN_ROOT", "RUN_ID", "OUTPUT_FILENAME"
    ))
    assert "UPSTREAM_PATH" not in params and "OUTPUT_PATH" not in joined
    assert "promote_run" not in final and "publish_compatibility_symlink" not in joined
    assert not any(
        isinstance(node, ast.Call) and getattr(node.func, "id", None) == "prepare_run"
        for source in code if source != final for node in ast.walk(ast.parse(source))
    )
    for marker in (
        "STABILITY_ENABLED", "adjusted_rand_score", "rank_genes_groups", "silhouette",
        "sc.pl.umap", "cluster_sizes", "doublet_score", "dendrogram", "pd.crosstab",
        "GASTRIC_MARKERS",
    ):
        assert marker in joined


@pytest.mark.parametrize(
    ("case", "match"),
    [
        ("ok", None), ("warnings", None), ("run_id", "run_id"), ("stage", "stage"),
        ("path", "escapes"), ("hash", "hash"),
        ("missing_status", "stage_status"), ("failed", "stage_status"),
        ("needs_review", "stage_status"), ("unaccepted_warnings", "warning_acceptance"),
    ],
)
def test_stage05_reads_only_verified_promoted_stage04(
    tmp_path: Path, case: str, match: str | None
) -> None:
    if match:
        with pytest.raises((ValueError, FileNotFoundError), match=match):
            _load_env(tmp_path, case)
    else:
        env, adata = _load_env(tmp_path, case)
        assert env["adata"] is adata
        upstream = env["upstream_input"]
        assert upstream["stage"] == "04_embedded"
        assert rc.sha256_file(upstream["manifest_path"]) == upstream["manifest_sha256"]
        assert rc.sha256_file(upstream["checkpoint_path"]) == upstream["checkpoint_sha256"]


def test_valid_stage05_is_auditable_unpromoted_draft(tmp_path: Path) -> None:
    env, adata = _final_env(tmp_path)
    exec(_cell("# === Stage 05 draft checkpoint"), env)
    manifest_path = tmp_path / "runs/stage05/draft/manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["stage_status"] == "NEEDS_REVIEW"
    assert manifest["candidate_cluster_keys"] == ["leiden_res_0.2", "leiden_res_0.4"]
    assert manifest["cluster_counts"] == {"leiden_res_0.2": 2, "leiden_res_0.4": 3}
    assert {"RESOLUTIONS", "USE_REP", "N_NEIGHBORS", "RUN_ID", "OUTPUT_FILENAME"} <= manifest["effective_parameters"].keys()
    assert manifest["effective_parameters"]["RESOLUTIONS"] == [0.2, 0.4]
    assert "UPSTREAM_CHECKPOINT" not in manifest["effective_parameters"]
    assert all(manifest["inputs"][0][key] for key in ("manifest_sha256", "checkpoint_sha256"))
    assert rc.validate_checkpoint(manifest_path).is_file()
    assert adata.uns["status"] == "NEEDS_REVIEW"
    assert not (tmp_path / "runs/stage05/promoted").exists()
    assert not (tmp_path / "results/05_clustered_v1.h5ad").exists()
    with pytest.raises(ValueError, match="cannot be promoted"):
        rc.promote_run(env["run_paths"])


def test_real_anndata_h5ad_roundtrip_preserves_run_metadata(tmp_path: Path) -> None:
    env, _ = _final_env(tmp_path)
    obs = pd.DataFrame({
        "leiden_res_0.2": ["0", "0", "1", "1"],
        "leiden_res_0.4": ["0", "1", "2", "2"],
    }, index=[f"cell-{i}" for i in range(4)])
    env["adata"] = ad.AnnData(sp.csr_matrix(np.ones((4, 2), dtype=np.float32)), obs=obs)
    exec(_cell("# === Stage 05 draft checkpoint"), env)
    restored = ad.read_h5ad(env["draft_output_path"])
    upstream = restored.uns["upstream_inputs"]["stage04"]
    assert (restored.uns["status"], restored.uns["run_id"]) == ("NEEDS_REVIEW", "stage05")
    assert upstream["stage"] == "04_embedded" and upstream["checkpoint_sha256"] == env["upstream_input"]["checkpoint_sha256"]
    assert rc.validate_checkpoint(env["run_paths"].manifest_path).is_file()


@pytest.mark.parametrize(
    ("case", "gate"),
    [
        ("missing", "candidate_columns_present"),
        ("nan", "candidate_columns_complete"),
        ("single", "candidate_columns_have_multiple_clusters"),
    ],
)
def test_invalid_candidate_columns_write_failed_audit(
    tmp_path: Path, case: str, gate: str
) -> None:
    env, adata = _final_env(tmp_path)
    if case == "missing":
        del adata.obs["leiden_res_0.4"]
    elif case == "nan":
        adata.obs.loc[0, "leiden_res_0.4"] = None
    else:
        adata.obs["leiden_res_0.4"] = "0"
    with pytest.raises(RuntimeError, match="Stage 05 FAILED"):
        exec(_cell("# === Stage 05 draft checkpoint"), env)
    draft = tmp_path / "runs/stage05/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["hard_postconditions"][gate] is False
    assert manifest["stage_status"] == "FAILED" and not list(draft.glob("*.h5ad"))


@pytest.mark.parametrize("filename", ["../x.h5ad", "manifest.json", "bad\nname.h5ad", "absolute"])
def test_unsafe_output_filename_never_writes_checkpoint(
    tmp_path: Path, filename: str
) -> None:
    env, adata = _final_env(tmp_path)
    if filename == "absolute":
        filename = str(tmp_path / "escaped.h5ad")
    env["OUTPUT_FILENAME"] = filename
    with pytest.raises(RuntimeError, match="Stage 05 FAILED"):
        exec(_cell("# === Stage 05 draft checkpoint"), env)
    draft = tmp_path / "runs/stage05/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["hard_postconditions"]["output_filename_valid"] is False
    assert adata.write_calls == 0 and list(draft.iterdir()) == [draft / "manifest.json"]


def test_write_failure_removes_partial_and_records_failure(tmp_path: Path) -> None:
    env, adata = _final_env(tmp_path)
    adata.fail_write = True
    with pytest.raises(OSError, match="disk full"):
        exec(_cell("# === Stage 05 draft checkpoint"), env)
    draft = tmp_path / "runs/stage05/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage_status"] == "FAILED"
    assert manifest["failure"] == {"type": "OSError", "message": "disk full"}
    assert "checkpoint" not in manifest and not list(draft.glob("*.h5ad"))


def _contract_env(tmp_path: Path) -> tuple[dict, _Adata]:
    """Set up env with adata ready for the USE_REP dimension + contract validation cell."""
    env, adata = _load_env(tmp_path)
    return env, adata


_EXPRESSION_CONTRACT_VALID = {
    "x_scale": "normalized_log1p",
    "counts_layer": "counts",
    "counts_source": ".raw.X",
    "counts_validated": True,
    "counts_integer_check": "full",
    "soupx_layer": None,
    "processing_history": [],
    "stage": "03",
}


def test_expression_contract_validated_at_stage05(tmp_path: Path) -> None:
    """Cell 5ed26e7e validates expression_contract with stage=03, normalized_log1p."""
    env, adata = _contract_env(tmp_path)
    adata.uns["expression_contract"] = dict(_EXPRESSION_CONTRACT_VALID)
    source = _cell("# === USE_REP 维度匹配检查 ===")
    exec(source, env)
    # No exception means validation passed


def test_expression_contract_missing_rejects(tmp_path: Path) -> None:
    """Missing expression_contract raises KeyError."""
    env, adata = _contract_env(tmp_path)
    source = _cell("# === USE_REP 维度匹配检查 ===")
    with pytest.raises((KeyError, NameError), match="expression_contract"):
        exec(source, env)


def test_expression_contract_wrong_scale_rejects(tmp_path: Path) -> None:
    """Non-normalized_log1p x_scale raises ValueError."""
    env, adata = _contract_env(tmp_path)
    adata.uns["expression_contract"] = {
        **_EXPRESSION_CONTRACT_VALID, "x_scale": "raw_counts",
    }
    source = _cell("# === USE_REP 维度匹配检查 ===")
    with pytest.raises(ValueError, match="x_scale"):
        exec(source, env)


def test_expression_contract_wrong_stage_rejects(tmp_path: Path) -> None:
    """Wrong stage value raises ValueError."""
    env, adata = _contract_env(tmp_path)
    adata.uns["expression_contract"] = {
        **_EXPRESSION_CONTRACT_VALID, "stage": "01",
    }
    source = _cell("# === USE_REP 维度匹配检查 ===")
    with pytest.raises(ValueError, match="stage"):
        exec(source, env)
