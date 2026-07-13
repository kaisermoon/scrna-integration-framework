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
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
os.chdir(PROJECT_ROOT)
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

PYTHON = sys.executable
JUPYTER = shutil.which("jupyter") or "jupyter"

# 执行产物输出目录：用系统临时目录（跨平台）而非硬编码 /tmp
OUTPUT_DIR = os.path.join(tempfile.gettempdir(), "smoke_run_notebooks_outputs")

# Notebook list in execution order
NOTEBOOKS = [
    # Stage 01 - per dataset (4)
    ("01_kim", "notebooks/01_per_dataset/01_kim.ipynb", "results/01_kim_v1.h5ad"),
    ("01_nancang", "notebooks/01_per_dataset/01_nancang.ipynb", "results/01_nancang_v1.h5ad"),
    ("01_nowicki", "notebooks/01_per_dataset/01_nowicki.ipynb", "results/01_nowicki_v1.h5ad"),
    ("01_yue", "notebooks/01_per_dataset/01_yue.ipynb", "results/01_yue_v1.h5ad"),
    # Stage 02-06 - core pipeline (5)
    ("02_merged", "notebooks/02_merged.ipynb", "results/02_merged_v1.h5ad"),
    ("03_normalized", "notebooks/03_normalized.ipynb", "results/03_normalized_v1.h5ad"),
    ("04_embedded", "notebooks/04_embedded.ipynb", "results/04_embedded_v1.h5ad"),
    ("05_clustered", "notebooks/05_clustered.ipynb", "results/05_clustered_v1.h5ad"),
    ("06_annotated", "notebooks/06_annotated.ipynb", "results/06_annotated_v1.h5ad"),
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


def run_notebook(name, rel_path, output_check):
    """Run a single notebook with nbconvert, return results."""
    nb_path = os.path.join(PROJECT_ROOT, rel_path)
    if not os.path.exists(nb_path):
        return {"name": name, "status": "MISSING", "error": f"NB not found: {nb_path}"}

    exec_path = nb_path.replace('.ipynb', f'_exec_{name}.ipynb')

    print(f"\n{'='*60}")
    print(f"Running: {name} ({rel_path})")
    print(f"{'='*60}")

    output_path = os.path.join(PROJECT_ROOT, output_check) if output_check else None
    output_before = None
    if output_path and os.path.exists(output_path):
        stat = os.stat(output_path)
        output_before = (stat.st_mtime_ns, stat.st_size)

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

    # Check if expected output file exists
    if output_path:
        output_exists = os.path.exists(output_path)
        result_info["output_exists"] = output_exists
        if output_exists:
            stat = os.stat(output_path)
            output_after = (stat.st_mtime_ns, stat.st_size)
            result_info["output_size_bytes"] = stat.st_size
            result_info["output_fresh"] = output_before is None or output_after != output_before
        else:
            result_info["output_fresh"] = False

        if notebook_status == "PASS" and not output_exists:
            notebook_status = result_info["status"] = "MISSING_OUTPUT"
        elif notebook_status == "PASS" and not result_info["output_fresh"]:
            notebook_status = result_info["status"] = "STALE_OUTPUT"

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
    for index, (name, rel_path, output_check) in enumerate(NOTEBOOKS):
        result = run_notebook(name, rel_path, output_check)
        all_results.append(result)
        if result["status"] != "PASS":
            core_failed = True
            blocked = NOTEBOOKS[index + 1:] + DOWNSTREAM
            all_results.extend(
                {
                    "name": blocked_name,
                    "status": "BLOCKED",
                    "time_seconds": 0,
                    "error": f"Required core stage {name} failed",
                }
                for blocked_name, _, _ in blocked
            )
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
            "output_exists": r.get("output_exists", None),
            "output_size_bytes": r.get("output_size_bytes", None),
            "returncode": r.get("returncode", None),
            "stderr_tail": r.get("stderr_tail", ""),
            "stdout_tail": r.get("stdout_tail", ""),
            "required": r.get("required", r["name"] in {n for n, _, _ in NOTEBOOKS}),
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
