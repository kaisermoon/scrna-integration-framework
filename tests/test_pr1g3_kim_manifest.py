from __future__ import annotations

import ast
import gc
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import yaml

import scrna_integration.run_contract as rc

ROOT = Path(__file__).parents[1]
NOTEBOOK = ROOT / "notebooks/01_per_dataset/01_kim.ipynb"


def _source(cell: dict) -> str:
    return "".join(cell["source"])


def _cell(marker: str) -> str:
    for cell in json.loads(NOTEBOOK.read_text(encoding="utf-8"))["cells"]:
        if marker in _source(cell):
            return _source(cell)
    raise AssertionError(marker)


def _write_inputs(root: Path, value: float = 1) -> Path:
    data, source = root / "data", root / "data/kim.h5ad"
    data.mkdir(parents=True, exist_ok=True)
    ad.AnnData(sp.csr_matrix(np.full((2, 2), value, dtype=np.float32)), obs=pd.DataFrame(index=["c1", "c2"]),
               var=pd.DataFrame(index=["GENE1", "GENE2"])).write_h5ad(source)
    (data / "manifest.yaml").write_text("source_dataset: Kim_2023\ninput:\n  format: h5ad\n  path: data/kim.h5ad\n", encoding="utf-8")
    return source


def _input_env(root: Path, run_id: str, reader=ad.read_h5ad) -> dict:
    env: dict = {}
    exec(_cell("# === PARAMS ==="), env)
    env.update(
        MANIFEST_PATH="data/manifest.yaml", RUN_ID=run_id, RUN_ROOT=str(root / "runs"),
        OUTPUT_FILENAME="kim.h5ad", _root=str(root), _started_at=rc.utc_now_rfc3339(),
        Path=Path, yaml=yaml, shutil=shutil, sp=sp, np=np, gc=gc,
        sc=SimpleNamespace(read_h5ad=reader), RSCRIPT_BIN="Rscript", R_AVAILABLE=False,
    )
    for name in (
        "MANIFEST_SCHEMA_VERSION", "MethodStatus", "artifact_record", "atomic_write_json",
        "collect_runtime_provenance", "determine_stage_status", "fingerprint_input",
        "prepare_run", "promote_run", "sha256_file", "snapshot_effective_parameters",
        "utc_now_rfc3339", "validate_manifest",
    ):
        env[name] = getattr(rc, name)
    exec(_cell("# === 数据读入"), env)
    return env


def _manifest(root: Path, run_id: str, state: str) -> tuple[Path, dict]:
    path = root / f"runs/{run_id}/{state}/manifest.json"
    return path, json.loads(path.read_text(encoding="utf-8"))


def test_kim_contract_cells_preserve_scientific_qc_and_managed_outputs() -> None:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    code = "\n".join(_source(cell) for cell in notebook["cells"] if cell["cell_type"] == "code")
    ast.parse(code)
    input_cell, final = _cell("# === 数据读入"), _cell("# Checkpoint")
    assert input_cell.index("_source_before") < input_cell.index("adata = sc.read_h5ad") < input_cell.index("_source_after") < input_cell.index("prepare_run")
    assert all(token in input_cell + final for token in ("fingerprint_input", "artifact_record", "validate_manifest", "schema_version"))
    assert "results/figures" not in code and code.count("_save_run_figure(") == 8
    assert all(_cell(marker) for marker in ("# 样本级 QC 摘要表", "# 自适应阈值计算", "# === N_MAD 敏感度分析", "# 过滤后 QC 小提琴图"))
    assert "adata.raw" not in input_cell + final and "layers[" not in input_cell + final


def test_actual_h5ad_bytes_change_input_sha_with_yaml_unchanged(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.chdir(tmp_path)
    source = _write_inputs(tmp_path)
    first = _input_env(tmp_path, "first")["_input_records"][1]["sha256"]
    yaml_before = (tmp_path / "data/manifest.yaml").read_bytes()
    _write_inputs(tmp_path, value=2)
    second = _input_env(tmp_path, "second")["_input_records"][1]["sha256"]
    assert source.is_file() and first != second and (tmp_path / "data/manifest.yaml").read_bytes() == yaml_before


@pytest.mark.parametrize("failure", ["mutation", "read"])
def test_input_failure_or_mutation_does_not_claim_run(monkeypatch, tmp_path: Path, failure: str) -> None:
    monkeypatch.chdir(tmp_path)
    _write_inputs(tmp_path)
    def reader(path):
        if failure == "read":
            raise OSError("read failed")
        value = ad.read_h5ad(path)
        Path(path).write_bytes(Path(path).read_bytes() + b"changed")
        return value

    with pytest.raises((OSError, RuntimeError)):
        _input_env(tmp_path, "rejected", reader)
    assert not (tmp_path / "runs").exists()


@pytest.mark.parametrize(("warnings", "review", "status", "state"), [
    ([], False, "SUCCESS", "promoted"), (["inspect QC"], False, "SUCCESS_WITH_WARNINGS", "draft"),
    ([], True, "NEEDS_REVIEW", "draft"),
])
def test_terminal_manifests_validate_with_live_artifacts(
    monkeypatch, tmp_path: Path, warnings, review, status: str, state: str
) -> None:
    monkeypatch.chdir(tmp_path)
    _write_inputs(tmp_path)
    env = _input_env(tmp_path, status.lower())
    env["STAGE_WARNINGS"], env["REQUIRE_REVIEW"] = warnings, review
    for name in ("01_kim_qc_violin_pre.png", "01_kim_qc_scatter_post.png"):
        (env["FIGURE_DIR"] / name).write_bytes(name.encode())
    exec(_cell("# Checkpoint"), env)
    path, manifest = _manifest(tmp_path, status.lower(), state)
    assert rc.validate_manifest(path)["stage_status"] == status
    checkpoint = rc.validate_checkpoint(path)
    assert set(rc.validate_artifacts(path)) == {item["role"] for item in manifest["artifacts"]}
    assert all(item["size"] == (path.parent / item["path"]).stat().st_size for item in manifest["artifacts"])
    saved = ad.read_h5ad(checkpoint)
    assert saved.uns["manifest_schema_version"] == rc.MANIFEST_SCHEMA_VERSION
    assert [item["role"] for item in manifest["inputs"]] == ["config_yaml", "source_h5ad"]
    assert all(set(item) == {"role", "path", "kind", "sha256"} and not Path(item["path"]).is_absolute()
               for item in manifest["inputs"])


@pytest.mark.parametrize("failure", ["gate", "write", "artifact"])
def test_failed_terminal_is_strict_and_cleans_partial(monkeypatch, tmp_path: Path, failure: str) -> None:
    monkeypatch.chdir(tmp_path)
    _write_inputs(tmp_path)
    env = _input_env(tmp_path, failure)
    (env["FIGURE_DIR"] / "diagnostic.png").write_bytes(b"diagnostic")
    if failure == "gate":
        env["adata"].X = np.ones((2, 2), dtype=np.float32)
    elif failure == "write":
        def fail_write(self, path, **kwargs):
            Path(path).write_bytes(b"partial")
            raise OSError("disk full")
        monkeypatch.setattr(ad.AnnData, "write_h5ad", fail_write)
    else:
        (env["FIGURE_DIR"] / "unsafe.png").symlink_to("diagnostic.png")
    with pytest.raises((RuntimeError, OSError, ValueError)):
        exec(_cell("# Checkpoint"), env)
    path, manifest = _manifest(tmp_path, failure, "draft")
    assert rc.validate_manifest(path)["stage_status"] == "FAILED"
    assert "checkpoint" not in manifest and not list(path.parent.glob("*.h5ad"))
    if failure == "artifact":
        assert manifest["artifacts"] == []
        assert list(env["FIGURE_DIR"].iterdir()) == []
    else:
        assert rc.validate_artifacts(path)["figure_diagnostic"].read_bytes() == b"diagnostic"


def test_artifact_directory_failure_writes_strict_failed_manifest(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.chdir(tmp_path)
    _write_inputs(tmp_path)
    mkdir = Path.mkdir
    def fail_figures(self, *args, **kwargs):
        if self.name == "figures":
            raise OSError("figure directory failed")
        return mkdir(self, *args, **kwargs)

    monkeypatch.setattr(Path, "mkdir", fail_figures)
    with pytest.raises(OSError, match="figure directory failed"):
        _input_env(tmp_path, "figure-failed")
    path, manifest = _manifest(tmp_path, "figure-failed", "draft")
    assert rc.validate_manifest(path)["stage_status"] == "FAILED"
    assert manifest["artifacts"] == [] and "checkpoint" not in manifest


def test_preterminal_qc_failure_leaves_auditable_sentinel(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.chdir(tmp_path)
    _write_inputs(tmp_path)
    env = _input_env(tmp_path, "qc-crash")
    figure = SimpleNamespace(savefig=lambda path, **kwargs: Path(path).write_bytes(b"qc"))
    env["_save_run_figure"](figure, "before_crash.png")
    env["adata"].var_names = ["GENE1", "gene2"]
    with pytest.raises(ValueError, match="基因 ID 轴不一致"):
        exec(_cell("# Checkpoint"), env)
    path, manifest = _manifest(tmp_path, "qc-crash", "draft")
    assert rc.validate_manifest(path)["failure"]["type"] == "IncompleteRun"
    assert rc.validate_artifacts(path)["figure_before_crash"].read_bytes() == b"qc"
    assert "checkpoint" not in manifest
