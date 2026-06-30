#!/usr/bin/env python3
"""Pipeline test runner: Execute each notebook, capture cell-level pass/fail, timing, errors."""

import shutil, subprocess, json, os, sys, time, re
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
os.chdir(PROJECT_ROOT)
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

PYTHON = sys.executable
JUPYTER = shutil.which("jupyter") or "jupyter"

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

# Downstream notebooks
DOWNSTREAM = [
    ("07_deg", "notebooks/07_downstream/07_deg.ipynb", ""),
    ("08_pseudobulk_deg", "notebooks/07_downstream/08_pseudobulk_deg.ipynb", ""),
    ("09_cnv", "notebooks/07_downstream/09_cnv.ipynb", ""),
    ("10_pseudotime", "notebooks/07_downstream/10_pseudotime.ipynb", ""),
    ("10b_pseudotime_monocle3", "notebooks/07_downstream/10b_pseudotime_monocle3.ipynb", ""),
    ("10c_pseudotime_cellrank2", "notebooks/07_downstream/10c_pseudotime_cellrank2.ipynb", ""),
    ("10d_pseudotime_cytotrace2", "notebooks/07_downstream/10d_pseudotime_cytotrace2.ipynb", ""),
    ("10e_pseudotime_compare", "notebooks/07_downstream/10e_pseudotime_compare.ipynb", ""),
    ("11_abundance", "notebooks/07_downstream/11_abundance.ipynb", ""),
    ("12_pathway", "notebooks/07_downstream/12_pathway.ipynb", ""),
    ("13_grn", "notebooks/07_downstream/13_grn.ipynb", ""),
    ("14_cell_communication", "notebooks/07_downstream/14_cell_communication.ipynb", ""),
    ("15_gene_modules", "notebooks/07_downstream/15_gene_modules.ipynb", ""),
    ("16_trajectory_de", "notebooks/07_downstream/16_trajectory_de.ipynb", ""),
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
            status = 'SKIP'  # No output = skipped/never executed
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
            status = f'FAIL'
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

    t_start = time.time()
    try:
        # Use absolute JUPYTER path to avoid picking up system jupyter-nbconvert
        # (PATH has /usr/local/bin before conda env, system one lacks pandocfilters)
        result = subprocess.run(
            [JUPYTER, "nbconvert",
             "--execute", "--to", "notebook",
             "--output", os.path.basename(exec_path),
             "--output-dir", "/tmp/pipeline_test_outputs",
             "--allow-errors",
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
    os.makedirs("/tmp/pipeline_test_outputs", exist_ok=True)
    nb_name = os.path.basename(nb_path).replace('.ipynb', f'_exec_{name}.ipynb')
    exec_out = os.path.join("/tmp/pipeline_test_outputs", nb_name)

    cells = []
    notebook_status = "ERROR"

    if os.path.exists(exec_out):
        cells = parse_cell_status(exec_out)
        if cells:
            n_pass = sum(1 for c in cells if c['status'] == 'PASS')
            n_fail = sum(1 for c in cells if c['status'] == 'FAIL')
            n_skip = sum(1 for c in cells if c['status'] == 'SKIP')
            if n_fail == 0:
                notebook_status = "PASS"
            else:
                notebook_status = "PARTIAL"
        else:
            notebook_status = "EMPTY"
    else:
        # Check stderr for clues
        notebook_status = "NO_OUTPUT"

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
    if output_check:
        out_path = os.path.join(PROJECT_ROOT, output_check)
        result_info["output_exists"] = os.path.exists(out_path)
        if os.path.exists(out_path):
            result_info["output_size_bytes"] = os.path.getsize(out_path)

    print(f"  Status: {notebook_status} | Time: {elapsed:.1f}s | Cells: {result_info['pass_cells']}P/{result_info['fail_cells']}F/{result_info['skip_cells']}S")
    return result_info


def main():
    os.makedirs("/tmp/pipeline_test_outputs", exist_ok=True)

    # Clean stale executed notebooks from previous runs
    out_dir = "/tmp/pipeline_test_outputs"
    for f in os.listdir(out_dir):
        if f.endswith("_executed.ipynb") or f.endswith("_exec.ipynb") or "_exec_" in f:
            os.remove(os.path.join(out_dir, f))

    all_results = []
    total_start = time.time()

    # Run core pipeline
    for name, rel_path, output_check in NOTEBOOKS:
        result = run_notebook(name, rel_path, output_check)
        all_results.append(result)

    # Run downstream notebooks (skip if core pipeline failed)
    for name, rel_path, output_check in DOWNSTREAM:
        result = run_notebook(name, rel_path, output_check)
        all_results.append(result)

    total_elapsed = time.time() - total_start

    # ===== Produce JSON summary =====
    summary = {
        "timestamp": datetime.now().isoformat(),
        "total_time_seconds": round(total_elapsed, 1),
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
        }
        # Add failed cell details
        if r.get("cells"):
            entry["failures"] = [
                {"cell_num": c["cell_num"], "error": c["error"], "source": c["source"]}
                for c in r["cells"] if c["status"] == "FAIL"
            ]
        summary["results"].append(entry)

    # Write summary JSON
    summary_path = "/tmp/pipeline_test_outputs/pipeline_test_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    # Print summary to stdout
    print(f"\n{'='*60}")
    print(f"PIPELINE TEST COMPLETE")
    print(f"Total time: {total_elapsed:.1f}s ({total_elapsed/60:.1f}min)")
    print(f"Summary written to: {summary_path}")
    print(f"{'='*60}")
    for r in all_results:
        status_icon = "OK" if r['status'] == 'PASS' else r['status']
        print(f"  {r['name']:25s}  {status_icon:10s}  {r.get('time_seconds', 0):8.1f}s  {r.get('pass_cells', 0)}P/{r.get('fail_cells', 0)}F/{r.get('skip_cells', 0)}S")


if __name__ == "__main__":
    main()
