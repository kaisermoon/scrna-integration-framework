#!/usr/bin/env python3
"""冒烟运行器：依次 nbconvert --execute 主线（01→06）与全部下游（D01→D14）notebook，
逐 cell 记录 pass/fail、耗时与报错，产出 JSON 汇总。

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

from scrna_integration.run_contract import sha256_file, validate_checkpoint

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
    # Stage 01 - per dataset (4)
    {"name": "01_kim", "notebook": "notebooks/01_per_dataset/01_kim.ipynb", "expected_stage": "01_qcd"},
    {"name": "01_nancang", "notebook": "notebooks/01_per_dataset/01_nancang.ipynb", "expected_stage": "01_qcd"},
    {"name": "01_nowicki", "notebook": "notebooks/01_per_dataset/01_nowicki.ipynb", "expected_stage": "01_qcd"},
    {"name": "01_yue", "notebook": "notebooks/01_per_dataset/01_yue.ipynb", "expected_stage": "01_qcd"},
    # Stage 02-06 - core pipeline (5)
    {"name": "02_merged", "notebook": "notebooks/02_merged.ipynb", "expected_stage": "02_merged"},
    {"name": "03_normalized", "notebook": "notebooks/03_normalized.ipynb", "expected_stage": "03_normalized"},
    {"name": "04_embedded", "notebook": "notebooks/04_embedded.ipynb", "expected_stage": "04_embedded"},
    {"name": "05_clustered", "notebook": "notebooks/05_clustered.ipynb", "expected_stage": "05_clustered"},
    {"name": "06_annotated", "notebook": "notebooks/06_annotated.ipynb", "expected_stage": "06_annotated"},
]

# Downstream notebooks（D 前缀命名：目录内各模块彼此无执行先后，均以
# 06_annotated.h5ad 为输入、可单独运行，D 编号仅为清单序位不表示顺序）
DOWNSTREAM = [
    ("D01_deg", "notebooks/07_downstream/D01_deg.ipynb", ""),
    ("D02_pseudobulk_deg", "notebooks/07_downstream/D02_pseudobulk_deg.ipynb", ""),
    ("D03_cnv", "notebooks/07_downstream/D03_cnv.ipynb", ""),
    ("D04_pseudotime", "notebooks/07_downstream/D04_pseudotime.ipynb", ""),
    ("D05_pseudotime_monocle3", "notebooks/07_downstream/D05_pseudotime_monocle3.ipynb", ""),
    ("D06_pseudotime_cellrank2", "notebooks/07_downstream/D06_pseudotime_cellrank2.ipynb", ""),
    ("D07_potency_cytotrace2", "notebooks/07_downstream/D07_potency_cytotrace2.ipynb", ""),
    ("D08_pseudotime_compare", "notebooks/07_downstream/D08_pseudotime_compare.ipynb", ""),
    ("D09_abundance", "notebooks/07_downstream/D09_abundance.ipynb", ""),
    ("D10_pathway", "notebooks/07_downstream/D10_pathway.ipynb", ""),
    ("D11_grn", "notebooks/07_downstream/D11_grn.ipynb", ""),
    ("D12_cell_communication", "notebooks/07_downstream/D12_cell_communication.ipynb", ""),
    ("D13_gene_modules", "notebooks/07_downstream/D13_gene_modules.ipynb", ""),
    ("D14_trajectory_de", "notebooks/07_downstream/D14_trajectory_de.ipynb", ""),
]


class ManifestValidationError(ValueError):
    """带稳定状态码的 runner manifest 验收错误。"""

    def __init__(self, code, message):
        self.code = code
        super().__init__(message)


_RUN_ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z", re.ASCII)


def snapshot_run_manifests(run_root=RUN_ROOT):
    """记录 run 根目录下 manifest 路径及内容 hash。"""
    root = Path(run_root)
    if not root.is_dir():
        return {}
    return {path: sha256_file(path) for path in root.glob("*/*/manifest.json") if path.is_file()}


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
    if status == "SUCCESS_WITH_WARNINGS":
        acceptance = manifest.get("warning_acceptance")
        if not isinstance(acceptance, dict) or not all(
            isinstance(acceptance.get(key), str) and acceptance[key].strip()
            for key in ("accepted_by", "accepted_at")
        ):
            raise ManifestValidationError("INVALID_MANIFEST", "warnings require accepted_by and accepted_at")
    if status == "NEEDS_REVIEW" and state != "draft":
        raise ManifestValidationError("INVALID_MANIFEST", "NEEDS_REVIEW must remain draft")
    if status != "NEEDS_REVIEW" and state != "promoted":
        raise ManifestValidationError("INVALID_MANIFEST", f"{status} must be promoted")
    try:
        checkpoint = validate_checkpoint(manifest_path)
    except (OSError, ValueError) as error:
        raise ManifestValidationError("INVALID_MANIFEST", f"checkpoint validation failed: {error}") from error
    _assert_manifest_hash(manifest_path, expected_hash)
    return {"run_id": run_id, "stage": manifest["stage"], "state": state, "status": status,
            "manifest": str(manifest_path), "checkpoint": str(checkpoint)}


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
        has_error = any(o.get('output_type') == 'error' for o in outputs)
        if not outputs:
            # Assignments and other valid cells often produce no display output.
            status = 'PASS' if cell.get('execution_count') is not None else 'SKIP'
        elif has_error:
            err_text = ''
            for o in outputs:
                if o.get('output_type') == 'error':
                    err_text = f"{o.get('ename', '')}: {o.get('evalue', '')}"
                    # Include traceback summary
                    traceback = o.get('traceback', [])
                    if traceback:
                        err_text += ' | ' + traceback[-1].strip()[:200]
                    break
            status = 'FAIL'
        else:
            status = 'PASS'
        cells.append({
            'cell_num': i + 1,
            'status': status,
            'error': err_text if has_error else '',
            'source': source
        })
    return cells


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
        return {"name": name, "status": "TIMEOUT", "time_seconds": elapsed, "error": "Execution timeout (>1h)"}

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

    if notebook_status == "PASS" and manifests_before is not None:
        try:
            result_info["run_outcome"] = validate_manifest_delta(
                manifests_before, RUN_ROOT, expected_stage
            )
        except ManifestValidationError as error:
            notebook_status = result_info["status"] = error.code
            result_info["error"] = str(error)

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

    # Core stages are required. Once one fails, downstream inputs are invalid.
    core_failed = False
    for index, notebook in enumerate(NOTEBOOKS):
        name = notebook["name"]
        result = run_notebook(name, notebook["notebook"], notebook["expected_stage"])
        all_results.append(result)
        if result["status"] != "PASS":
            core_failed = True
            blocked_names = [item["name"] for item in NOTEBOOKS[index + 1:]]
            blocked_names.extend(item[0] for item in DOWNSTREAM)
            all_results.extend({"name": blocked_name, "status": "BLOCKED", "time_seconds": 0,
                                "error": f"Required core stage {name} failed"}
                               for blocked_name in blocked_names)
            break

    # Downstream notebooks are optional diagnostics in this PR. Their failures
    # remain visible in the report but do not override a successful core run.
    if not core_failed:
        for name, rel_path, output_check in DOWNSTREAM:
            result = run_notebook(name, rel_path, output_check)
            result["required"] = False
            all_results.append(result)

    total_elapsed = time.time() - total_start
    optional_failed = any(
        r.get("required") is False and r["status"] != "PASS" for r in all_results
    )
    if core_failed:
        overall_status = "FAILED"
    elif optional_failed:
        overall_status = "SUCCESS_WITH_WARNINGS"
    else:
        overall_status = "SUCCESS"

    # ===== Produce JSON summary =====
    summary = {
        "timestamp": datetime.now().isoformat(),
        "total_time_seconds": round(total_elapsed, 1),
        "status": overall_status,
        "results": []
    }

    for r in all_results:
        entry = {
            "name": r["name"],
            "status": r["status"],
            "time_seconds": r.get("time_seconds", 0),
            "cells_total": r.get("total_cells", 0),
            "cells_pass": r.get("pass_cells", 0),
            "cells_fail": r.get("fail_cells", 0),
            "cells_skip": r.get("skip_cells", 0),
            "returncode": r.get("returncode", None),
            "stderr_tail": r.get("stderr_tail", ""),
            "stdout_tail": r.get("stdout_tail", ""),
            "required": r.get("required", r["name"] in {item["name"] for item in NOTEBOOKS}),
            "run_outcome": r.get("run_outcome"),
        }
        # Add failed cell details
        if r.get("cells"):
            entry["failures"] = [
                {"cell_num": c["cell_num"], "error": c["error"], "source": c["source"]}
                for c in r["cells"] if c["status"] == "FAIL"
            ]
        summary["results"].append(entry)

    # Write summary JSON
    summary_path = os.path.join(OUTPUT_DIR, "smoke_run_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    # Print summary to stdout
    print(f"\n{'='*60}")
    status_message = {
        "FAILED": "PIPELINE TEST FAILED",
        "SUCCESS_WITH_WARNINGS": "PIPELINE TEST COMPLETE WITH WARNINGS",
        "SUCCESS": "PIPELINE TEST COMPLETE",
    }[overall_status]
    print(status_message)
    print(f"Total time: {total_elapsed:.1f}s ({total_elapsed/60:.1f}min)")
    print(f"Summary written to: {summary_path}")
    print(f"{'='*60}")
    for r in all_results:
        status_icon = "OK" if r['status'] == 'PASS' else r['status']
        print(f"  {r['name']:25s}  {status_icon:10s}  {r.get('time_seconds', 0):8.1f}s  {r.get('pass_cells', 0)}P/{r.get('fail_cells', 0)}F/{r.get('skip_cells', 0)}S")

    return 1 if core_failed else 0


if __name__ == "__main__":
    sys.exit(main())
