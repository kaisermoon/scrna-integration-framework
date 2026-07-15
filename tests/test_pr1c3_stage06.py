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
NOTEBOOK = ROOT / "notebooks/06_annotated.ipynb"


def _source(cell: dict) -> str:
    return "".join(value) if isinstance(value := cell.get("source", ""), list) else value


def _cell(marker: str) -> str:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return next(_source(cell) for cell in notebook["cells"]
                if cell["cell_type"] == "code" and marker in _source(cell))


class _Adata:
    def __init__(self) -> None:
        self.X = sp.csr_matrix(np.ones((4, 2), dtype=np.float32))
        self.n_obs, self.n_vars = self.X.shape
        self.obs = pd.DataFrame({"leiden_res_0.6": ["0", "0", "1", "1"]})
        self.obsm, self.uns = {}, {}
        self.fail_write = False
        self.write_calls = 0

    def write_h5ad(self, path: Path, **_: object) -> None:
        self.write_calls += 1
        Path(path).write_bytes(b"partial" if self.fail_write else b"stage06")
        if self.fail_write:
            raise OSError("disk full")


def _upstream(root: Path) -> tuple[Path, Path]:
    paths = rc.prepare_run(root, "stage05")
    checkpoint = paths.draft_dir / "05.h5ad"
    checkpoint.write_bytes(b"stage05")
    rc.atomic_write_json(paths.manifest_path, {
        "run_id": "stage05", "stage": "05_clustered", "stage_status": "SUCCESS",
        "checkpoint": {"path": checkpoint.name, "sha256": rc.sha256_file(checkpoint)},
    })
    rc.promote_run(paths)
    return paths.promoted_dir / checkpoint.name, paths.promoted_dir / "manifest.json"


def _load_env(tmp_path: Path, case: str = "ok", cluster_case: str = "ok") -> tuple[dict, _Adata]:
    checkpoint, manifest_path = _upstream(tmp_path / "upstream")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if case == "run_id":
        manifest["run_id"] = "other"
    elif case == "stage":
        manifest["stage"] = "04_embedded"
    elif case == "path":
        manifest["checkpoint"]["path"] = "../05.h5ad"
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
    if cluster_case == "nonempty":
        adata.X, adata.n_obs, adata.n_vars, adata.obs = sp.csr_matrix((0, 0), dtype=np.float32), 0, 0, adata.obs.iloc[:0]
    elif cluster_case == "missing":
        del adata.obs["leiden_res_0.6"]
    elif cluster_case == "nan":
        adata.obs.loc[0, "leiden_res_0.6"] = None
    elif cluster_case == "empty":
        adata.obs["leiden_res_0.6"] = ""
    env = {
        "Path": Path, "json": json, "sp": sp, "atomic_write_json": rc.atomic_write_json,
        "prepare_run": rc.prepare_run, "resume_run": rc.resume_run,
        "sha256_file": rc.sha256_file, "validate_checkpoint": rc.validate_checkpoint,
        "snapshot_effective_parameters": rc.snapshot_effective_parameters,
        "collect_runtime_provenance": lambda *_: {"python": "test"},
        "sc": SimpleNamespace(read_h5ad=lambda path: adata if Path(path) == checkpoint else None),
    }
    exec(_cell("# === PARAMS ==="), env)
    env.update(UPSTREAM_RUN_ROOT=str(tmp_path / "upstream"), UPSTREAM_RUN_ID="stage05",
               RUN_ROOT=str(tmp_path / "runs"), RUN_ID="stage06", _root=str(tmp_path))
    setup = _cell("# === setup")
    start = setup.index("# === 只读取 verified promoted Stage 05 上游 ===")
    end = setup.index("# === Stage 06 input preflight 完成 ===")
    exec(setup[start:end], env)
    return env, adata


def _final_env(tmp_path: Path) -> tuple[dict, _Adata]:
    env, adata = _load_env(tmp_path)
    env.update(
        Path=Path, re=re, RUN_ROOT=str(tmp_path / "runs"), RUN_ID="stage06",
        OUTPUT_FILENAME="06.h5ad", _root=str(tmp_path),
        atomic_write_json=rc.atomic_write_json, collect_runtime_provenance=lambda *_: {"python": "test"},
        determine_stage_status=rc.determine_stage_status, prepare_run=rc.prepare_run,
        snapshot_effective_parameters=rc.snapshot_effective_parameters,
    )
    return env, adata


def test_structure_preserves_annotation_workflow_and_defers_pr6() -> None:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    code = [_source(cell) for cell in notebook["cells"] if cell["cell_type"] == "code"]
    joined = "\n".join(code)
    ast.parse(joined, filename=str(NOTEBOOK))
    params = _cell("# === PARAMS ===")
    setup = _cell("# === setup")
    final = _cell("# === Stage 06 draft checkpoint")
    assert all(name in params for name in ("UPSTREAM_RUN_ROOT", "UPSTREAM_RUN_ID", "RUN_ROOT", "RUN_ID", "OUTPUT_FILENAME"))
    assert "UPSTREAM_PATH" not in params and "OUTPUT_PATH" not in joined
    assert "promote_run" not in final and "PI_CONFIRMED" in joined
    assert "publish_compatibility_symlink" not in joined
    assert "if not all(input_hard_postconditions.values())" in setup
    assert not any(
        isinstance(node, ast.Call) and getattr(node.func, "id", None) == "prepare_run"
        for source in code if source not in {setup, final} for node in ast.walk(ast.parse(source))
    )
    for marker in ("marker_assignments", "mLLMCelltype", "score_genes", "scANVI",
                   "_verdict_results", "PI 拍板", "pi_decisions"):
        assert marker in joined


@pytest.mark.parametrize(("case", "match"), [
    ("ok", None), ("warnings", None), ("run_id", "run_id"), ("stage", "stage"),
    ("path", "escapes"), ("hash", "hash"), ("missing_status", "stage_status"),
    ("failed", "stage_status"), ("needs_review", "stage_status"),
    ("unaccepted_warnings", "warning_acceptance"),
])
def test_stage06_reads_only_verified_promoted_stage05(tmp_path: Path, case: str, match: str | None) -> None:
    if match:
        with pytest.raises((ValueError, FileNotFoundError), match=match):
            _load_env(tmp_path, case)
    else:
        env, adata = _load_env(tmp_path, case)
        upstream = env["upstream_input"]
        assert env["adata"] is adata and upstream["stage"] == "05_clustered"
        assert rc.sha256_file(upstream["manifest_path"]) == upstream["manifest_sha256"]
        assert rc.sha256_file(upstream["checkpoint_path"]) == upstream["checkpoint_sha256"]
        assert not (tmp_path / "runs/stage06").exists()


@pytest.mark.parametrize(("case", "gate"), [
    ("nonempty", "non_empty"), ("missing", "cluster_column_present"), ("nan", "cluster_labels_complete"),
    ("empty", "cluster_labels_nonempty"),
])
def test_input_preflight_audits_before_downstream_consumers(tmp_path: Path, case: str, gate: str) -> None:
    with pytest.raises(RuntimeError, match="input preflight FAILED"):
        _load_env(tmp_path, cluster_case=case)
    draft = tmp_path / "runs/stage06/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["hard_postconditions"][gate] is False
    assert manifest["stage_status"] == "FAILED" and manifest["failure"]["type"] == "InputContractError"
    assert not list(draft.glob("*.h5ad"))


def test_valid_stage06_is_auditable_unpromoted_draft(tmp_path: Path) -> None:
    env, adata = _final_env(tmp_path)
    exec(_cell("# === Stage 06 draft checkpoint"), env)
    manifest_path = tmp_path / "runs/stage06/draft/manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["stage_status"] == "NEEDS_REVIEW" and manifest["cluster_key"] == "leiden_res_0.6"
    assert {"LEIDEN_COL", "ANNOTATION_OUTPUT_VERSION", "RUN_ID", "OUTPUT_FILENAME"} <= manifest["effective_parameters"].keys()
    assert "UPSTREAM_CHECKPOINT" not in manifest["effective_parameters"]
    assert all(manifest["inputs"][0][key] for key in ("manifest_sha256", "checkpoint_sha256"))
    assert rc.validate_checkpoint(manifest_path).is_file()
    assert adata.uns["status"] == "NEEDS_REVIEW" and not (tmp_path / "runs/stage06/promoted").exists()
    with pytest.raises(ValueError, match="cannot be promoted"):
        rc.promote_run(env["run_paths"])


def test_real_anndata_h5ad_roundtrip_preserves_run_metadata(tmp_path: Path) -> None:
    env, _ = _final_env(tmp_path)
    obs = pd.DataFrame({"leiden_res_0.6": ["0", "0", "1", "1"]},
                       index=[f"cell-{i}" for i in range(4)])
    env["adata"] = ad.AnnData(sp.csr_matrix(np.ones((4, 2), dtype=np.float32)), obs=obs)
    exec(_cell("# === Stage 06 draft checkpoint"), env)
    restored = ad.read_h5ad(env["draft_output_path"])
    upstream = restored.uns["upstream_inputs"]["stage05"]
    assert (restored.uns["status"], restored.uns["run_id"]) == ("NEEDS_REVIEW", "stage06")
    assert upstream["stage"] == "05_clustered" and upstream["checkpoint_sha256"] == env["upstream_input"]["checkpoint_sha256"]
    assert rc.validate_checkpoint(env["run_paths"].manifest_path).is_file()


@pytest.mark.parametrize(("case", "gate"), [
    ("nonempty", "non_empty"),
    ("missing", "cluster_column_present"),
    ("nan", "cluster_labels_complete"),
    ("empty", "cluster_labels_nonempty"),
])
def test_invalid_cluster_labels_write_failed_audit(
    tmp_path: Path, case: str, gate: str
) -> None:
    env, adata = _final_env(tmp_path)
    if case == "nonempty":
        adata.n_obs = adata.n_vars = 0
    elif case == "missing":
        del adata.obs["leiden_res_0.6"]
    elif case == "nan":
        adata.obs.loc[0, "leiden_res_0.6"] = None
    else:
        adata.obs["leiden_res_0.6"] = ""
    with pytest.raises(RuntimeError, match="Stage 06 FAILED"):
        exec(_cell("# === Stage 06 draft checkpoint"), env)
    draft = tmp_path / "runs/stage06/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["hard_postconditions"][gate] is False
    assert manifest["stage_status"] == "FAILED" and not list(draft.glob("*.h5ad"))


@pytest.mark.parametrize("filename", ["../x.h5ad", "manifest.json", "bad\nname.h5ad", "absolute"])
def test_unsafe_output_filename_never_writes_checkpoint(tmp_path: Path, filename: str) -> None:
    env, adata = _final_env(tmp_path)
    if filename == "absolute":
        filename = str(tmp_path / "escaped.h5ad")
    env["OUTPUT_FILENAME"] = filename
    with pytest.raises(RuntimeError, match="Stage 06 FAILED"):
        exec(_cell("# === Stage 06 draft checkpoint"), env)
    draft = tmp_path / "runs/stage06/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["hard_postconditions"]["output_filename_valid"] is False
    assert adata.write_calls == 0 and list(draft.iterdir()) == [draft / "manifest.json"]


def test_write_failure_cleans_partial_and_run_collision_is_refused(tmp_path: Path) -> None:
    env, adata = _final_env(tmp_path)
    adata.fail_write = True
    with pytest.raises(OSError, match="disk full"):
        exec(_cell("# === Stage 06 draft checkpoint"), env)
    draft = tmp_path / "runs/stage06/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["failure"] == {"type": "OSError", "message": "disk full"}
    assert manifest["stage_status"] == "FAILED" and not list(draft.glob("*.h5ad"))
    collision_env, _ = _final_env(tmp_path / "collision")
    rc.prepare_run(collision_env["RUN_ROOT"], collision_env["RUN_ID"])
    with pytest.raises(FileExistsError):
        exec(_cell("# === Stage 06 draft checkpoint"), collision_env)
