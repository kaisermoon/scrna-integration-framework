"""PR 0: smoke/CI failures must remain observable to callers."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = ROOT / "scripts" / "smoke_run_notebooks.py"


def _load_script(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _load_runner():
    return _load_script("smoke_run_notebooks", RUNNER_PATH)


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

    result = runner.run_notebook("broken", "broken.ipynb")

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

    result = runner.run_notebook("broken", "broken.ipynb")

    assert result["status"] == "ERROR"
    assert result["returncode"] == 2


def test_runner_optional_run_ignores_missing_legacy_output(monkeypatch, tmp_path: Path) -> None:
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

    result = runner.run_notebook("stage", "stage.ipynb")

    assert result["status"] == "PASS"
    assert "output_exists" not in result


def test_runner_optional_run_ignores_unchanged_legacy_output(monkeypatch, tmp_path: Path) -> None:
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

    result = runner.run_notebook("stage", "stage.ipynb")

    assert result["status"] == "PASS"
    assert "output_fresh" not in result


def test_core_failure_blocks_downstream_and_returns_nonzero(monkeypatch, tmp_path: Path) -> None:
    runner = _load_runner()
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path))
    monkeypatch.setattr(runner, "NOTEBOOKS", [
        {"name": "core", "notebook": "core.ipynb", "expected_stage": "core"},
        {"name": "later", "notebook": "later.ipynb", "expected_stage": "later"},
    ])
    calls = []

    def fake_run(name, rel_path, output_check):
        calls.append(name)
        return {"name": name, "status": "ERROR", "time_seconds": 0}

    monkeypatch.setattr(runner, "run_notebook", fake_run)

    assert runner.main() == 1
    assert calls == ["core"]
    summary = json.loads((tmp_path / "smoke_run_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "FAILED"
    assert summary["results"][1]["status"] == "BLOCKED_FAILURE"


def test_action_required_stops_core_and_returns_two(monkeypatch, tmp_path: Path) -> None:
    runner = _load_runner()
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path))
    monkeypatch.setattr(runner, "NOTEBOOKS", [
        {"name": "04_embedded", "notebook": "04.ipynb", "expected_stage": "04_embedded"},
        {"name": "05_clustered", "notebook": "05.ipynb", "expected_stage": "05_clustered"},
    ])
    calls = []

    def fake_run(name, rel_path, output_check):
        calls.append(name)
        return {"name": name, "stage": "04_embedded", "status": "REVIEW_REQUIRED", "time_seconds": 0}

    monkeypatch.setattr(runner, "run_notebook", fake_run)

    assert runner.main() == 2
    assert calls == ["04_embedded"]
    summary = json.loads((tmp_path / "smoke_run_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "ACTION_REQUIRED"
    assert summary["results"][1]["status"] == "BLOCKED_REVIEW"


def test_all_pass_has_success_overall_status(monkeypatch, tmp_path: Path) -> None:
    runner = _load_runner()
    monkeypatch.setattr(runner, "OUTPUT_DIR", str(tmp_path))
    monkeypatch.setattr(runner, "NOTEBOOKS", [
        {"name": "core", "notebook": "core.ipynb", "expected_stage": "core"}])
    monkeypatch.setattr(
        runner,
        "run_notebook",
        lambda name, rel_path, output_check: {"name": name, "status": "PASS"},
    )

    assert runner.main() == 0
    summary = json.loads((tmp_path / "smoke_run_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "SUCCESS"


def test_r_check_builds_parse_command_and_propagates_failure(monkeypatch, tmp_path: Path) -> None:
    r_check = _load_script("check_r_scripts", ROOT / "scripts" / "check_r_scripts.py")
    r_file = tmp_path / "entry point.R"
    calls = []
    monkeypatch.setattr(
        r_check.subprocess,
        "run",
        lambda command, **kwargs: calls.append((command, kwargs)) or SimpleNamespace(returncode=7),
    )
    monkeypatch.setattr(
        sys, "argv", ["check_r_scripts.py", "--rscript", "custom-Rscript", str(r_file)]
    )

    assert r_check.main() == 7
    parse_expr = "for (f in commandArgs(trailingOnly=TRUE)) parse(file=f)"
    expected = ["custom-Rscript", "--vanilla", "-e", parse_expr, str(r_file)]
    assert calls == [(expected, {"check": False})]
