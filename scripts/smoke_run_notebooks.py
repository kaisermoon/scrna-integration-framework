#!/usr/bin/env python3
"""冒烟运行器：依次执行并验收主线（01→06）notebook，
逐 cell 记录 pass/fail、耗时、人工动作点与报错，产出 JSON 汇总。

这不是 pytest 用例（不被 tests/ 收集），是端到端手动冒烟脚本，用于验证
notebook 之间的 h5ad 输入输出契约不断裂。命令行直接运行：

    python scripts/smoke_run_notebooks.py

执行产物（各 notebook 的已执行副本 + 汇总 JSON）写入系统临时目录下的
smoke_run_notebooks_outputs/，不污染仓库。
"""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

from scrna_integration.run_contract import atomic_write_json, sha256_file, validate_checkpoint

PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
os.chdir(PROJECT_ROOT)
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

PYTHON = sys.executable
JUPYTER = shutil.which("jupyter") or "jupyter"

# 执行产物输出目录：用系统临时目录（跨平台）而非硬编码 /tmp
OUTPUT_DIR = os.path.join(tempfile.gettempdir(), "smoke_run_notebooks_outputs")
RUN_ROOT = Path(PROJECT_ROOT) / "results" / "runs"

# Notebook list in execution order
NOTEBOOKS = [
    # Stage 01 - 按输入数据格式命名的四个模板（4）
    # 01 系列按「输入数据格式」而非「数据集」命名：数据集是下游项目的资产，
    # 格式才是框架该提供的可复用骨架。模板默认指向 data/_subset/ 夹具。
    {"name": "01_template_10x_mtx", "notebook": "notebooks/01_per_dataset/01_template_10x_mtx.ipynb", "expected_stage": "01_qcd"},
    {"name": "01_template_10x_h5", "notebook": "notebooks/01_per_dataset/01_template_10x_h5.ipynb", "expected_stage": "01_qcd"},
    {"name": "01_template_h5ad", "notebook": "notebooks/01_per_dataset/01_template_h5ad.ipynb", "expected_stage": "01_qcd"},
    {"name": "01_template_counts_matrix", "notebook": "notebooks/01_per_dataset/01_template_counts_matrix.ipynb", "expected_stage": "01_qcd"},
    # Stage 02-06 - core pipeline (5)
    {"name": "02_merged", "notebook": "notebooks/02_merged.ipynb", "expected_stage": "02_merged"},
    {"name": "03_normalized", "notebook": "notebooks/03_normalized.ipynb", "expected_stage": "03_normalized"},
    {"name": "04_embedded", "notebook": "notebooks/04_embedded.ipynb", "expected_stage": "04_embedded"},
    {"name": "05_clustered", "notebook": "notebooks/05_clustered.ipynb", "expected_stage": "05_clustered"},
    {"name": "06_annotated", "notebook": "notebooks/06_annotated.ipynb", "expected_stage": "06_annotated"},
]

class ManifestValidationError(ValueError):
    """带稳定状态码的 runner manifest 验收错误。"""

    def __init__(self, code, message):
        self.code = code
        super().__init__(message)


_RUN_ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z", re.ASCII)
MAX_FAILED_CELLS = 20
MAX_ENAME_CHARS, MAX_EVALUE_CHARS, MAX_TRACEBACK_CHARS = 100, 500, 1000
TRUNCATION_MARKER = "...[truncated]"


def _bounded_error(value, limit):
    text = str(value).strip()
    return text if len(text) <= limit else text[:limit - len(TRUNCATION_MARKER)] + TRUNCATION_MARKER


def snapshot_run_manifests(run_root=RUN_ROOT):
    """记录 run 根目录下 manifest 路径及内容 hash。"""
    root = Path(run_root)
    if not root.is_dir():
        return {}
    return {path: sha256_file(path) for path in root.glob("*/*/manifest.json") if path.is_file()}


def snapshot_run_dirs(run_root=RUN_ROOT):
    root = Path(run_root)
    return set(root.iterdir()) if root.is_dir() else set()


def _assert_manifest_hash(manifest_path, expected_hash):
    try:
        current_hash = sha256_file(manifest_path)
    except OSError as error:
        raise ManifestValidationError("STALE_RUN", f"new manifest disappeared: {error}") from error
    if current_hash != expected_hash:
        raise ManifestValidationError("STALE_RUN", "new manifest changed during validation")


def _manifest_location(manifest_path, run_root):
    try:
        relative = Path(manifest_path).relative_to(Path(run_root))
    except ValueError as error:
        raise ManifestValidationError("INVALID_MANIFEST", "manifest is outside run root") from error
    if len(relative.parts) != 3 or relative.name != "manifest.json":
        raise ManifestValidationError(
            "INVALID_MANIFEST", "manifest path must be <run_id>/<draft|promoted>/manifest.json"
        )
    run_id, state, _ = relative.parts
    if _RUN_ID_RE.fullmatch(run_id) is None:
        raise ManifestValidationError("INVALID_MANIFEST", f"unsafe run_id directory: {run_id!r}")
    if state not in {"draft", "promoted"}:
        raise ManifestValidationError("INVALID_MANIFEST", f"unknown run state: {state}")
    return run_id, state


def validate_manifest_delta(before, run_root=RUN_ROOT, expected_stage=None):
    """验收一次 core 执行唯一新增且完整的 run manifest。"""
    after = snapshot_run_manifests(run_root)
    changed = [path for path in before if after.get(path) != before[path]]
    if changed:
        names = ", ".join(str(path) for path in sorted(changed))
        raise ManifestValidationError("STALE_RUN", f"preexisting manifest was modified: {names}")
    new_paths = sorted(after.keys() - before.keys())
    if not new_paths:
        raise ManifestValidationError("MISSING_MANIFEST", "execution produced no new run manifest")
    if len(new_paths) != 1:
        raise ManifestValidationError(
            "MULTIPLE_MANIFESTS", f"execution produced {len(new_paths)} new run manifests"
        )
    manifest_path = new_paths[0]
    expected_hash = after[manifest_path]
    _assert_manifest_hash(manifest_path, expected_hash)
    run_id, state = _manifest_location(manifest_path, run_root)
    try:
        with manifest_path.open(encoding="utf-8") as stream:
            manifest = json.load(stream)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ManifestValidationError("INVALID_MANIFEST", f"cannot read manifest: {error}") from error
    if not isinstance(manifest, dict):
        raise ManifestValidationError("INVALID_MANIFEST", "manifest must be a JSON object")
    if manifest.get("run_id") != run_id:
        raise ManifestValidationError("INVALID_MANIFEST", "manifest run_id does not match directory")
    if manifest.get("stage") != expected_stage:
        raise ManifestValidationError(
            "INVALID_MANIFEST", f"expected stage {expected_stage}, got {manifest.get('stage')}"
        )
    status = manifest.get("stage_status")
    if status not in {"SUCCESS", "SUCCESS_WITH_WARNINGS", "NEEDS_REVIEW"}:
        raise ManifestValidationError("INVALID_MANIFEST", f"unacceptable stage status: {status}")
    acceptance_valid = False
    if status == "SUCCESS_WITH_WARNINGS":
        acceptance = manifest.get("warning_acceptance")
        acceptance_valid = isinstance(acceptance, dict) and all(
            isinstance(acceptance.get(key), str) and acceptance[key].strip()
            for key in ("accepted_by", "accepted_at")
        )
        if state == "promoted" and not acceptance_valid:
            raise ManifestValidationError("INVALID_MANIFEST", "warnings require accepted_by and accepted_at")
    if status == "NEEDS_REVIEW" and state != "draft":
        raise ManifestValidationError("INVALID_MANIFEST", "NEEDS_REVIEW must remain draft")
    try:
        checkpoint = validate_checkpoint(manifest_path)
    except (OSError, ValueError) as error:
        raise ManifestValidationError("INVALID_MANIFEST", f"checkpoint validation failed: {error}") from error
    _assert_manifest_hash(manifest_path, expected_hash)
    if state == "promoted":
        action = "PASS"
    elif status == "NEEDS_REVIEW" or (status == "SUCCESS_WITH_WARNINGS" and not acceptance_valid):
        action = "REVIEW_REQUIRED"
    else:
        action = "READY_TO_PROMOTE"
    return {"run_id": run_id, "stage": manifest["stage"], "state": state, "status": status,
            "action": action, "manifest": str(manifest_path), "checkpoint": str(checkpoint)}


def audit_partial_run(before_dirs, expected_stage, failure, run_root=RUN_ROOT):
    """为唯一安全 draft 写审计，再尽力清理；审计错误不覆盖原始失败。"""
    try:
        return _audit_partial_run(before_dirs, expected_stage, failure, run_root)
    except Exception as error:
        return {"run_id": None, "manifest": None,
                "audit_error": f"unexpected audit failure: {error}"}


def _audit_partial_run(before_dirs, expected_stage, failure, run_root):
    if before_dirs is None or Path(run_root).is_symlink():
        return None
    new_dirs = [path for path in snapshot_run_dirs(run_root) - before_dirs
                if path.is_dir() or path.is_symlink()]
    if len(new_dirs) != 1:
        return None
    run_dir = new_dirs[0]
    if run_dir.is_symlink() or _RUN_ID_RE.fullmatch(run_dir.name) is None:
        return None
    draft = run_dir / "draft"
    manifest_path = draft / "manifest.json"
    if not draft.is_dir() or draft.is_symlink() or manifest_path.exists() or manifest_path.is_symlink():
        return None
    if (run_dir / "promoted").exists() or (run_dir / "promoted").is_symlink():
        return None
    payload = {"run_id": run_dir.name, "stage": expected_stage,
               "stage_status": "FAILED", "failure": failure}
    try:
        atomic_write_json(manifest_path, payload)
    except Exception as error:
        return {"run_id": run_dir.name, "manifest": None,
                "audit_error": f"audit write failed: {error}"}
    cleanup_errors = []
    for child in list(draft.iterdir()):
        if child == manifest_path:
            continue
        try:
            if child.is_symlink() or child.is_file():
                child.unlink()
            elif child.is_dir():
                shutil.rmtree(child)
            else:
                raise ValueError("unsupported partial entry")
        except Exception as error:
            cleanup_errors.append(f"{child.name}: {error}")
    audit_error = None
    if cleanup_errors:
        payload["audit_cleanup_errors"] = cleanup_errors
        audit_error = "; ".join(cleanup_errors)
        try:
            atomic_write_json(manifest_path, payload, overwrite=True)
        except Exception as error:
            audit_error += f"; cleanup error persistence failed: {error}"
    return {"run_id": run_dir.name, "manifest": str(manifest_path),
            "audit_error": audit_error}


def parse_cell_status(executed_nb_path):
    """Parse an executed notebook JSON for per-cell status."""
    with open(executed_nb_path) as f:
        nb = json.load(f)
    cells = []
    for i, cell in enumerate(nb['cells']):
        if cell['cell_type'] != 'code':
            continue
        source = ''.join(cell['source'])[:120].replace('\n', ' ')
        outputs = cell.get('outputs', [])
        error_output = next((o for o in outputs if o.get('output_type') == 'error'), None)
        has_error = error_output is not None
        ename = evalue = traceback_summary = ''
        if not outputs:
            # Assignments and other valid cells often produce no display output.
            status = 'PASS' if cell.get('execution_count') is not None else 'SKIP'
        elif has_error:
            ename = _bounded_error(error_output.get('ename', ''), MAX_ENAME_CHARS)
            evalue = _bounded_error(error_output.get('evalue', ''), MAX_EVALUE_CHARS)
            traceback = error_output.get('traceback', [])
            traceback_summary = _bounded_error(traceback[-1], MAX_TRACEBACK_CHARS) if traceback else ''
            err_text = f"{ename}: {evalue}" + (f" | {traceback_summary}" if traceback_summary else '')
            status = 'FAIL'
        else:
            status = 'PASS'
        cells.append({
            'cell_num': i + 1,
            'status': status,
            'error': err_text if has_error else '',
            'ename': ename,
            'evalue': evalue,
            'traceback_summary': traceback_summary,
            'source': source
        })
    return cells


def _failed_cell_records(cells):
    failed = [cell for cell in cells if cell['status'] == 'FAIL']
    records = [{"cell_index": cell["cell_num"], "ename": cell["ename"], "evalue": cell["evalue"],
                "traceback_summary": cell["traceback_summary"]} for cell in failed[:MAX_FAILED_CELLS]]
    return records, max(0, len(failed) - MAX_FAILED_CELLS)


def run_notebook(name, rel_path, expected_stage=None):
    """Run a single notebook with nbconvert, return results."""
    nb_path = os.path.join(PROJECT_ROOT, rel_path)
    if not os.path.exists(nb_path):
        return {"name": name, "status": "MISSING", "error": f"NB not found: {nb_path}"}

    exec_path = nb_path.replace('.ipynb', f'_exec_{name}.ipynb')

    print(f"\n{'='*60}")
    print(f"Running: {name} ({rel_path})")
    print(f"{'='*60}")

    manifests_before = snapshot_run_manifests(RUN_ROOT) if expected_stage else None
    run_dirs_before = snapshot_run_dirs(RUN_ROOT) if expected_stage else None

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    t_start = time.time()
    try:
        # Use absolute JUPYTER path to avoid picking up system jupyter-nbconvert
        # (PATH has /usr/local/bin before conda env, system one lacks pandocfilters)
        result = subprocess.run(
            [JUPYTER, "nbconvert",
             "--execute", "--to", "notebook",
             "--output", os.path.basename(exec_path),
             "--output-dir", OUTPUT_DIR,
             "--ExecutePreprocessor.timeout=3600",
             nb_path],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=3700,
            env={**os.environ, "PYTHONPATH": f"{PROJECT_ROOT}/src:{os.environ.get('PYTHONPATH', '')}"}
        )
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t_start
        audit = audit_partial_run(run_dirs_before, expected_stage,
                                  {"returncode": None, "stderr": "Execution timeout (>1h)"}, RUN_ROOT)
        return {"name": name, "stage": expected_stage, "status": "TIMEOUT",
                "time_seconds": elapsed, "error": "Execution timeout (>1h)",
                "run_id": audit.get("run_id") if audit else None,
                "manifest": audit.get("manifest") if audit else None,
                "audit_error": audit.get("audit_error") if audit else None}

    elapsed = time.time() - t_start

    # Check if executed notebook was produced
    nb_name = os.path.basename(nb_path).replace('.ipynb', f'_exec_{name}.ipynb')
    exec_out = os.path.join(OUTPUT_DIR, nb_name)

    cells = []
    notebook_status = "ERROR"

    if os.path.exists(exec_out):
        cells = parse_cell_status(exec_out)
        if cells:
            n_fail = sum(1 for c in cells if c['status'] == 'FAIL')
            if result.returncode != 0 or n_fail:
                notebook_status = "ERROR"
            else:
                notebook_status = "PASS"
        else:
            notebook_status = "EMPTY"
    else:
        notebook_status = "ERROR" if result.returncode != 0 else "NO_OUTPUT"

    if result.returncode != 0:
        notebook_status = "ERROR"

    result_info = {
        "name": name,
        "stage": expected_stage,
        "status": notebook_status,
        "time_seconds": round(elapsed, 1),
        "total_cells": len(cells),
        "pass_cells": sum(1 for c in cells if c['status'] == 'PASS'),
        "fail_cells": sum(1 for c in cells if c['status'] == 'FAIL'),
        "skip_cells": sum(1 for c in cells if c['status'] == 'SKIP'),
        "stdout_tail": result.stdout[-500:] if result.stdout else "",
        "stderr_tail": result.stderr[-500:] if result.stderr else "",
        "returncode": result.returncode,
        "cells": cells,
        "executed_nb_path": exec_out
    }
    if notebook_status != "PASS":
        failed_cell = next((cell for cell in cells if cell['status'] == 'FAIL'), None)
        result_info["error"] = (failed_cell["error"] if failed_cell else
                                result_info["stderr_tail"] or f"Notebook execution status: {notebook_status}")

    if notebook_status == "PASS" and manifests_before is not None:
        try:
            outcome = validate_manifest_delta(manifests_before, RUN_ROOT, expected_stage)
            result_info["run_outcome"] = outcome
            notebook_status = result_info["status"] = result_info["action"] = outcome["action"]
        except ManifestValidationError as error:
            notebook_status = result_info["status"] = error.code
            result_info["error"] = str(error)
    elif notebook_status != "PASS" and run_dirs_before is not None:
        failure = {"returncode": result.returncode}
        if result.stderr:
            failure["stderr"] = result.stderr[-500:]
        if os.path.exists(exec_out):
            failure["executed_notebook"] = exec_out
        cell_errors, omitted = _failed_cell_records(cells)
        if cell_errors:
            failure["cell_errors"] = cell_errors
            failure["truncated_cell_error_count"] = omitted
        audit = audit_partial_run(run_dirs_before, expected_stage, failure, RUN_ROOT)
        if audit:
            result_info.update(audit)

    print(f"  Status: {notebook_status} | Time: {elapsed:.1f}s | Cells: {result_info['pass_cells']}P/{result_info['fail_cells']}F/{result_info['skip_cells']}S")
    return result_info


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Clean stale executed notebooks from previous runs
    out_dir = OUTPUT_DIR
    for f in os.listdir(out_dir):
        if f.endswith("_executed.ipynb") or f.endswith("_exec.ipynb") or "_exec_" in f:
            os.remove(os.path.join(out_dir, f))

    all_results = []
    total_start = time.time()

    exit_code = 0
    overall_status = "SUCCESS"
    for index, notebook in enumerate(NOTEBOOKS):
        name = notebook["name"]
        result = run_notebook(name, notebook["notebook"], notebook["expected_stage"])
        all_results.append(result)
        if result["status"] != "PASS":
            action_required = result["status"] in {"READY_TO_PROMOTE", "REVIEW_REQUIRED"}
            exit_code = 2 if action_required else 1
            overall_status = "ACTION_REQUIRED" if action_required else "FAILED"
            blocked_status = "BLOCKED_REVIEW" if action_required else "BLOCKED_FAILURE"
            all_results.extend({"name": item["name"], "stage": item["expected_stage"],
                                "status": blocked_status, "action": blocked_status,
                                "time_seconds": 0, "error": f"Core stage {name} stopped the run"}
                               for item in NOTEBOOKS[index + 1:])
            break

    total_elapsed = time.time() - total_start

    # ===== Produce JSON summary =====
    summary = {
        "timestamp": datetime.now().isoformat(),
        "total_time_seconds": round(total_elapsed, 1),
        "status": overall_status,
        "results": []
    }

    for r in all_results:
        outcome = r.get("run_outcome") or {}
        entry = {
            "name": r["name"],
            "stage": outcome.get("stage", r.get("stage")),
            "run_id": outcome.get("run_id", r.get("run_id")),
            "manifest": outcome.get("manifest", r.get("manifest")),
            "status": r["status"],
            "action": r.get("action", outcome.get("action")),
            "error": r.get("error"),
            "audit_error": r.get("audit_error"),
            "time_seconds": r.get("time_seconds", 0),
            "cells_total": r.get("total_cells", 0),
            "cells_pass": r.get("pass_cells", 0),
            "cells_fail": r.get("fail_cells", 0),
            "cells_skip": r.get("skip_cells", 0),
            "returncode": r.get("returncode", None),
            "stderr_tail": r.get("stderr_tail", ""),
            "stdout_tail": r.get("stdout_tail", ""),
        }
        # Add failed cell details
        if r.get("cells"):
            entry["failures"], entry["truncated_cell_error_count"] = _failed_cell_records(r["cells"])
        summary["results"].append(entry)

    # Write summary JSON
    summary_path = os.path.join(OUTPUT_DIR, "smoke_run_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    # Print summary to stdout
    print(f"\n{'='*60}")
    status_message = {
        "FAILED": "PIPELINE TEST FAILED",
        "ACTION_REQUIRED": "PIPELINE ACTION REQUIRED",
        "SUCCESS": "PIPELINE TEST COMPLETE",
    }[overall_status]
    print(status_message)
    print(f"Total time: {total_elapsed:.1f}s ({total_elapsed/60:.1f}min)")
    print(f"Summary written to: {summary_path}")
    print(f"{'='*60}")
    for r in all_results:
        status_icon = "OK" if r['status'] == 'PASS' else r['status']
        print(f"  {r['name']:25s}  {status_icon:10s}  {r.get('time_seconds', 0):8.1f}s  {r.get('pass_cells', 0)}P/{r.get('fail_cells', 0)}F/{r.get('skip_cells', 0)}S")

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
