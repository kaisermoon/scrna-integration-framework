import ast
import gc
import json
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import scipy.sparse as sp

import scrna_integration.run_contract as rc

ROOT = Path(__file__).parents[1]
NOTEBOOK = ROOT / "notebooks/04_embedded.ipynb"


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


class _Adata:
    def __init__(self, embedding: np.ndarray | None = None) -> None:
        self.X = sp.csr_matrix(np.ones((2, 3), dtype=np.float32))
        self.n_obs, self.n_vars = self.X.shape
        self.obsm = {} if embedding is None else {"X_pca": embedding}
        self.uns: dict = {}
        self.layers: dict = {}
        self.fail_write = False
        self.write_calls = 0

    def write_h5ad(self, path: Path, **_: object) -> None:
        self.write_calls += 1
        Path(path).write_bytes(b"partial" if self.fail_write else b"checkpoint")
        if self.fail_write:
            raise OSError("disk full")


def _promoted_stage03(root: Path, adata: _Adata) -> tuple[Path, str]:
    run_id = "stage03"
    paths = rc.prepare_run(root, run_id)
    checkpoint = paths.draft_dir / "03.h5ad"
    checkpoint.write_bytes(b"upstream")
    rc.atomic_write_json(
        paths.manifest_path,
        {
            "run_id": run_id,
            "stage": "03_normalized",
            "stage_status": "SUCCESS",
            "checkpoint": {"path": checkpoint.name, "sha256": rc.sha256_file(checkpoint)},
        },
    )
    promoted = rc.promote_run(paths)
    return promoted, run_id


def _env(tmp_path: Path, embedding: np.ndarray | None) -> tuple[dict, _Adata]:
    adata = _Adata(embedding)
    upstream_root = tmp_path / "upstream"
    checkpoint, upstream_id = _promoted_stage03(upstream_root, adata)
    env = {
        "Path": Path,
        "json": json,
        "np": np,
        "re": re,
        "sp": sp,
        "gc": gc,
        "resume_run": rc.resume_run,
        "sha256_file": rc.sha256_file,
        "validate_checkpoint": rc.validate_checkpoint,
        "sc": SimpleNamespace(
            read_h5ad=lambda path: adata if Path(path) == checkpoint else None
        ),
        "UPSTREAM_RUN_ROOT": str(upstream_root),
        "UPSTREAM_RUN_ID": upstream_id,
    }
    setup = _cell("# === 只读取 promoted Stage 03 上游 ===")
    start = setup.index("# === 只读取 promoted Stage 03 上游 ===")
    end = setup.index("# === promoted Stage 03 上游读取完成 ===")
    exec(setup[start:end], env)
    exec(_cell("# === PARAMS ==="), env)
    env.update(
        RUN_ROOT=str(tmp_path / "runs"),
        RUN_ID="stage04",
        OUTPUT_FILENAME="04.h5ad",
        _root=str(tmp_path),
        atomic_write_json=rc.atomic_write_json,
        collect_runtime_provenance=lambda _, packages: {
            "packages_requested": sorted(packages)
        },
        determine_stage_status=rc.determine_stage_status,
        prepare_run=rc.prepare_run,
        snapshot_effective_parameters=rc.snapshot_effective_parameters,
        upstream_input=env["upstream_input"],
    )
    return env, adata


def test_notebook_contract_is_deferred_and_does_not_select_or_promote() -> None:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    code = [_source(cell) for cell in notebook["cells"] if cell["cell_type"] == "code"]
    ast.parse("\n".join(code), filename=str(NOTEBOOK))
    params = _cell("# === PARAMS ===")
    final = _cell("# === Stage 04 draft checkpoint")
    assert all(name in params for name in ("RUN_ROOT", "RUN_ID", "OUTPUT_FILENAME"))
    assert "OUTPUT_PATH" not in params and "SELECTED_EMBEDDING" not in "\n".join(code)
    assert "promote_run" not in final
    assert not any(
        isinstance(node, ast.Call) and getattr(node.func, "id", None) == "prepare_run"
        for source in code
        if source != final
        for node in ast.walk(ast.parse(source))
    )


def test_valid_stage04_is_auditable_unpromoted_draft(tmp_path: Path) -> None:
    env, adata = _env(tmp_path, np.ones((2, 3), dtype=np.float32))
    exec(_cell("# === Stage 04 draft checkpoint"), env)
    manifest_path = tmp_path / "runs/stage04/draft/manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["stage_status"] == "NEEDS_REVIEW"
    assert manifest["inputs"][0]["stage"] == "03_normalized"
    assert manifest["embedding_shapes"] == {"X_pca": [2, 3]}
    assert set(manifest["runtime_provenance"]["packages_requested"]) == {
        "anndata", "scanpy", "numpy", "pandas", "scipy", "harmonypy",
        "scvi-tools", "torch", "scikit-learn", "sccraft",
    }
    assert "BATCH_KEY" in manifest["effective_parameters"]
    assert "UPSTREAM_CHECKPOINT" not in manifest["effective_parameters"]
    assert rc.validate_checkpoint(manifest_path).is_file()
    assert adata.uns["status"] == "NEEDS_REVIEW"
    assert not (tmp_path / "runs/stage04/promoted").exists()
    with pytest.raises(FileNotFoundError):
        rc.resume_run(tmp_path / "runs", "stage04", promoted=True)
    with pytest.raises(ValueError, match="cannot be promoted"):
        rc.promote_run(env["run_paths"])


@pytest.mark.parametrize(
    ("embedding", "failed_gate"),
    [
        (None, "embedding_present"),
        (np.ones((3, 2), dtype=np.float32), "embedding_rows_match"),
        (np.empty((2, 0), dtype=np.float32), "embedding_dimensions_valid"),
        (np.array([[1.0, np.nan], [2.0, 3.0]], dtype=np.float32), "embeddings_finite"),
    ],
)
def test_invalid_embedding_writes_failed_audit_without_checkpoint(
    tmp_path: Path, embedding: np.ndarray | None, failed_gate: str
) -> None:
    env, _ = _env(tmp_path, embedding)
    with pytest.raises(RuntimeError, match="Stage 04 FAILED"):
        exec(_cell("# === Stage 04 draft checkpoint"), env)
    draft = tmp_path / "runs/stage04/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage_status"] == "FAILED"
    assert manifest["hard_postconditions"][failed_gate] is False
    assert "checkpoint" not in manifest and not list(draft.glob("*.h5ad"))


@pytest.mark.parametrize(
    "output_filename", ["../x.h5ad", "absolute", "manifest.json", "bad\nname.h5ad"]
)
def test_invalid_output_filename_cannot_write_outside_draft(
    tmp_path: Path, output_filename: str
) -> None:
    env, adata = _env(tmp_path, np.ones((2, 3), dtype=np.float32))
    if output_filename == "absolute":
        output_filename = str(tmp_path / "escaped.h5ad")
    env["OUTPUT_FILENAME"] = output_filename
    draft = tmp_path / "runs/stage04/draft"
    escaped = (draft / output_filename).resolve()
    with pytest.raises(RuntimeError, match="Stage 04 FAILED"):
        exec(_cell("# === Stage 04 draft checkpoint"), env)
    manifest_path = draft / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["hard_postconditions"]["output_filename_valid"] is False
    assert manifest["stage_status"] == "FAILED" and "checkpoint" not in manifest
    assert adata.write_calls == 0
    if escaped != manifest_path:
        assert not escaped.exists()
    assert list(draft.iterdir()) == [manifest_path]


def test_write_failure_removes_partial_checkpoint_and_records_failure(tmp_path: Path) -> None:
    env, adata = _env(tmp_path, np.ones((2, 3), dtype=np.float32))
    adata.fail_write = True
    with pytest.raises(OSError, match="disk full"):
        exec(_cell("# === Stage 04 draft checkpoint"), env)
    draft = tmp_path / "runs/stage04/draft"
    manifest = json.loads((draft / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage_status"] == "FAILED"
    assert manifest["failure"] == {"type": "OSError", "message": "disk full"}
    assert "checkpoint" not in manifest and not list(draft.glob("*.h5ad"))


def test_custom_embedding_is_covered_by_finiteness_gate(tmp_path: Path) -> None:
    env, adata = _env(tmp_path, np.ones((2, 3), dtype=np.float32))
    adata.obsm = {"X_custom": np.array([[1.0], [np.inf]], dtype=np.float32)}
    with pytest.raises(RuntimeError, match="Stage 04 FAILED"):
        exec(_cell("# === Stage 04 draft checkpoint"), env)
    manifest = json.loads(
        (tmp_path / "runs/stage04/draft/manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["embedding_shapes"] == {"X_custom": [2, 1]}
    assert manifest["hard_postconditions"]["embeddings_finite"] is False
