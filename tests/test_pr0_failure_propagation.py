"""PR 0: smoke/CI failures must remain observable to callers."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = ROOT / "scripts" / "smoke_run_notebooks.py"
R_CHECK_PATH = ROOT / "scripts" / "check_r_scripts.py"


def _load_runner():
    spec = importlib.util.spec_from_file_location("smoke_run_notebooks", RUNNER_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _load_r_check():
    spec = importlib.util.spec_from_file_location("check_r_scripts", R_CHECK_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _write_executed_notebook(path: Path, *, error: bool = False) -> None:
    outputs = []
    if error:
        outputs = [{"output_type": "error", "ename": "RuntimeError", "evalue": "boom"}]
    path.write_text(
        json.dumps(
            {
                "cells": [
                    {
                        "cell_type": "code",
                        "source": ["raise RuntimeError('boom')" if error else "x = 1"],
                        "execution_count": 1,
                        "outputs": outputs,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )


def test_runner_rejects_python_cell_error(monkeypatch, tmp_path: Path) -> None:
    runner = _load_runner()
    notebook = tmp_path / "broken.ipynb"
    notebook.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(runner, "PROJECT_ROOT", str(tmp_path))
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path / "executed"))

    def fake_run(*args, **kwargs):
        output = Path(runner.OUTPUT_DIR) / "broken_exec_broken.ipynb"
        output.parent.mkdir(exist_ok=True)
        _write_executed_notebook(output, error=True)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(runner.subprocess, "run", fake_run)

    result = runner.run_notebook("broken", "broken.ipynb", "")

    assert result["status"] == "ERROR"
    assert result["fail_cells"] == 1


def test_runner_rejects_nonzero_nbconvert(monkeypatch, tmp_path: Path) -> None:
    runner = _load_runner()
    notebook = tmp_path / "broken.ipynb"
    notebook.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(runner, "PROJECT_ROOT", str(tmp_path))
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path / "executed"))
    monkeypatch.setattr(
        runner.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=2, stdout="", stderr="failed"),
    )

    result = runner.run_notebook("broken", "broken.ipynb", "")

    assert result["status"] == "ERROR"
    assert result["returncode"] == 2


def test_runner_rejects_missing_required_output(monkeypatch, tmp_path: Path) -> None:
    runner = _load_runner()
    notebook = tmp_path / "stage.ipynb"
    notebook.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(runner, "PROJECT_ROOT", str(tmp_path))
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path / "executed"))

    def fake_run(*args, **kwargs):
        output = Path(runner.OUTPUT_DIR) / "stage_exec_stage.ipynb"
        output.parent.mkdir(exist_ok=True)
        _write_executed_notebook(output)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(runner.subprocess, "run", fake_run)

    result = runner.run_notebook("stage", "stage.ipynb", "results/required.h5ad")

    assert result["status"] == "MISSING_OUTPUT"
    assert result["output_exists"] is False


def test_runner_rejects_unchanged_required_output(monkeypatch, tmp_path: Path) -> None:
    runner = _load_runner()
    notebook = tmp_path / "stage.ipynb"
    notebook.write_text("{}", encoding="utf-8")
    required = tmp_path / "results" / "required.h5ad"
    required.parent.mkdir()
    required.write_text("stale", encoding="utf-8")
    monkeypatch.setattr(runner, "PROJECT_ROOT", str(tmp_path))
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path / "executed"))

    def fake_run(*args, **kwargs):
        output = Path(runner.OUTPUT_DIR) / "stage_exec_stage.ipynb"
        output.parent.mkdir(exist_ok=True)
        _write_executed_notebook(output)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(runner.subprocess, "run", fake_run)

    result = runner.run_notebook("stage", "stage.ipynb", "results/required.h5ad")

    assert result["status"] == "STALE_OUTPUT"
    assert result["output_fresh"] is False


def test_core_failure_blocks_downstream_and_returns_nonzero(monkeypatch, tmp_path: Path) -> None:
    runner = _load_runner()
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path))
    monkeypatch.setattr(runner, "NOTEBOOKS", [("core", "core.ipynb", "required.h5ad")])
    monkeypatch.setattr(runner, "DOWNSTREAM", [("optional", "optional.ipynb", "")])
    calls = []

    def fake_run(name, rel_path, output_check):
        calls.append(name)
        return {"name": name, "status": "ERROR", "time_seconds": 0}

    monkeypatch.setattr(runner, "run_notebook", fake_run)

    assert runner.main() == 1
    assert calls == ["core"]
    summary = json.loads((tmp_path / "smoke_run_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "FAILED"
    assert summary["results"][1]["status"] == "BLOCKED"


def test_optional_failure_is_reported_without_blocking_core(monkeypatch, tmp_path: Path) -> None:
    runner = _load_runner()
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path))
    monkeypatch.setattr(runner, "NOTEBOOKS", [("core", "core.ipynb", "required.h5ad")])
    monkeypatch.setattr(runner, "DOWNSTREAM", [("optional", "optional.ipynb", "")])

    def fake_run(name, rel_path, output_check):
        status = "PASS" if name == "core" else "ERROR"
        return {"name": name, "status": status, "time_seconds": 0}

    monkeypatch.setattr(runner, "run_notebook", fake_run)

    assert runner.main() == 0
    summary = json.loads((tmp_path / "smoke_run_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "SUCCESS_WITH_WARNINGS"
    assert summary["results"][1]["status"] == "ERROR"
    assert summary["results"][1]["required"] is False


def test_all_pass_has_success_overall_status(monkeypatch, tmp_path: Path) -> None:
    runner = _load_runner()
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path))
    monkeypatch.setattr(runner, "NOTEBOOKS", [("core", "core.ipynb", "required.h5ad")])
    monkeypatch.setattr(runner, "DOWNSTREAM", [("optional", "optional.ipynb", "")])
    monkeypatch.setattr(
        runner,
        "run_notebook",
        lambda name, rel_path, output_check: {"name": name, "status": "PASS"},
    )

    assert runner.main() == 0
    summary = json.loads((tmp_path / "smoke_run_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "SUCCESS"


def test_r_check_builds_parse_command(monkeypatch, tmp_path: Path) -> None:
    r_check = _load_r_check()
    rscript = tmp_path / "custom-Rscript"
    r_file = tmp_path / "entry point.R"
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=7)

    monkeypatch.setattr(r_check.subprocess, "run", fake_run)
    monkeypatch.setattr(
        sys,
        "argv",
        ["check_r_scripts.py", "--rscript", str(rscript), str(r_file)],
    )

    assert r_check.main() == 7
    assert calls == [
        (
            [
                str(rscript),
                "--vanilla",
                "-e",
                "for (f in commandArgs(trailingOnly=TRUE)) parse(file=f)",
                str(r_file),
            ],
            {"check": False},
        )
    ]


def test_r_check_propagates_rscript_failure(tmp_path: Path) -> None:
    fake_rscript = tmp_path / "Rscript"
    fake_rscript.write_text("#!/bin/sh\nexit 7\n", encoding="utf-8")
    fake_rscript.chmod(0o755)
    r_file = tmp_path / "broken.R"
    r_file.write_text("stop('boom')\n", encoding="utf-8")

    result = subprocess.run(
        [sys.executable, str(R_CHECK_PATH), "--rscript", str(fake_rscript), str(r_file)],
        check=False,
    )

    assert result.returncode == 7
