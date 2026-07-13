import ast
import gc
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from scrna_integration.run_contract import (
    atomic_write_json,
    determine_stage_status,
    prepare_run,
    promote_run,
    resume_run,
    sha256_file,
    validate_checkpoint,
)


PROJECT_ROOT = Path(__file__).parents[1]
NOTEBOOKS = [
    "notebooks/01_per_dataset/01_kim.ipynb",
    "notebooks/01_per_dataset/01_nancang.ipynb",
    "notebooks/01_per_dataset/01_nowicki.ipynb",
    "notebooks/01_per_dataset/01_yue.ipynb",
    "notebooks/02_merged.ipynb",
    "notebooks/03_normalized.ipynb",
]


def _notebook(path: str) -> dict:
    return json.loads((PROJECT_ROOT / path).read_text(encoding="utf-8"))


def _cell(path: str, marker: str) -> str:
    for cell in _notebook(path)["cells"]:
        source = "".join(cell.get("source", []))
        if marker in source:
            return source
    raise AssertionError(f"{path} 缺少 cell: {marker}")


class _FakeAdata:
    def __init__(self, n_obs: int = 2) -> None:
        self.X = sp.csr_matrix(np.ones((n_obs, 2), dtype=np.float32))
        self.n_obs, self.n_vars = self.X.shape
        self.var_names = ["GENE1", "GENE2"]
        self.obs = pd.DataFrame({"source_dataset": ["test"] * n_obs})
        self.var = pd.DataFrame({"highly_variable": [True, False]})
        self.uns: dict = {}

    def write_h5ad(self, path: Path, **_: object) -> None:
        Path(path).write_bytes(b"test-checkpoint")


def _stage01_env(tmp_path: Path, path: str, run_id: str, *, n_obs: int = 2) -> dict:
    env: dict = {}
    exec(_cell(path, "# === PARAMS ==="), env)
    manifest = tmp_path / "manifest.yaml"
    manifest.write_text("source_dataset: test\n", encoding="utf-8")
    env.update(
        adata=_FakeAdata(n_obs),
        MANIFEST_PATH=str(manifest),
        RUN_ID=run_id,
        RUN_ROOT=str(tmp_path / "runs"),
        run_paths=prepare_run(tmp_path / "runs", run_id),
        os=os,
        sp=sp,
        np=np,
        Path=Path,
        atomic_write_json=atomic_write_json,
        determine_stage_status=determine_stage_status,
        promote_run=promote_run,
        sha256_file=sha256_file,
        gc=gc,
    )
    return env


def test_all_01_03_notebooks_parse_and_expose_run_contract() -> None:
    for path in NOTEBOOKS:
        code = "\n".join(
            "".join(cell.get("source", []))
            for cell in _notebook(path)["cells"]
            if cell["cell_type"] == "code"
        )
        ast.parse(code, filename=path)
        assert "RUN_ID" in code and "RUN_ROOT" in code
        assert "prepare_run" in code and "promote_run" in code
        assert "stage_status" in code and "sha256" in code

    for path in NOTEBOOKS[-2:]:
        code = "\n".join("".join(c.get("source", [])) for c in _notebook(path)["cells"])
        assert "UPSTREAM_RUN_ID" in code
        assert "validate_checkpoint" in code


@pytest.mark.parametrize(
    ("path", "run_id"), [(path, f"stage01-source-{i}") for i, path in enumerate(NOTEBOOKS[:4])]
)
def test_stage01_checkpoint_promotes_success_and_run_id_cannot_overwrite(
    tmp_path: Path, path: str, run_id: str
) -> None:
    env = _stage01_env(tmp_path, path, run_id)
    exec(_cell(path, "# Checkpoint"), env)
    promoted = tmp_path / "runs" / run_id / "promoted"
    manifest = json.loads((promoted / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage_status"] == "SUCCESS"
    assert manifest["inputs"][0]["sha256"] == sha256_file(env["MANIFEST_PATH"])
    assert manifest["effective_parameters"]["qc_strategy"] == env["QC_STRATEGY"]
    assert validate_checkpoint(promoted / "manifest.json").is_file()
    with pytest.raises(FileExistsError):
        prepare_run(tmp_path / "runs", run_id)


def test_stage01_failed_gate_never_promotes(tmp_path: Path) -> None:
    env = _stage01_env(tmp_path, NOTEBOOKS[0], "stage01-failed", n_obs=0)
    with pytest.raises(RuntimeError, match="FAILED"):
        exec(_cell(NOTEBOOKS[0], "# Checkpoint：写入 per-dataset h5ad"), env)
    assert not (tmp_path / "runs" / "stage01-failed" / "promoted").exists()


def test_stage02_rejects_tampered_promoted_upstream(tmp_path: Path) -> None:
    upstream = prepare_run(tmp_path / "upstream", "stage01-source")
    checkpoint = upstream.draft_dir / "01_source.h5ad"
    checkpoint.write_bytes(b"valid")
    atomic_write_json(
        upstream.manifest_path,
        {
            "run_id": upstream.run_id,
            "stage_status": "SUCCESS",
            "checkpoint": {"path": checkpoint.name, "sha256": sha256_file(checkpoint)},
        },
    )
    promoted = promote_run(upstream)
    promoted.write_bytes(b"tampered")
    env = {
        "UPSTREAM_RUN_ROOT": str(tmp_path / "upstream"),
        "UPSTREAM_RUN_IDS": ["stage01-source"],
        "resume_run": resume_run,
        "validate_checkpoint": validate_checkpoint,
    }
    with pytest.raises(ValueError, match="hash"):
        exec(_cell(NOTEBOOKS[4], "# 加载每个 per-dataset h5ad"), env)
