#!/usr/bin/env python
"""Run 07_downstream/09_cnv.ipynb programmatically."""

import nbformat, sys, os, traceback, time
from pathlib import Path

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
os.chdir(_PROJECT_ROOT)

# Ensure src/ is on path (notebook does this too)
_src_path = os.path.join(os.getcwd(), 'src')
if _src_path not in sys.path:
    sys.path.insert(0, _src_path)

os.makedirs('results/figures', exist_ok=True)
os.makedirs('results/tables', exist_ok=True)

nb_path = 'notebooks/07_downstream/09_cnv.ipynb'
# Load notebook
nb = nbformat.read(nb_path, as_version=4)
code_cells = [(i, c) for i, c in enumerate(nb.cells) if c.cell_type == 'code']

print(f'=== {nb_path} ({len(code_cells)} code cells) ===')

results = []
for idx, (cell_idx, cell) in enumerate(code_cells):
    cell_num = idx + 1
    src = cell.source
    # Skip empty cells
    if not src.strip():
        results.append((cell_num, cell_idx, 'SKIP', 'empty'))
        continue

    print(f'[{cell_num}/{len(code_cells)}] cell[{cell_idx}]...', end=' ', flush=True)
    t0 = time.time()
    try:
        exec(src)
        elapsed = time.time() - t0
        print(f'OK ({elapsed:.1f}s)')
        results.append((cell_num, cell_idx, 'PASS', f'{elapsed:.1f}s'))
    except Exception as e:
        elapsed = time.time() - t0
        err = traceback.format_exc()
        print(f'FAIL ({elapsed:.1f}s)')
        # Print last few lines of traceback
        lines = err.strip().split('\n')
        for l in lines[-5:]:
            print(f'  {l}')
        results.append((cell_num, cell_idx, 'FAIL', str(e)[:200]))
        # Continue to next cell (don't stop)

# Summary
passed = sum(1 for r in results if r[2] == 'PASS')
failed = sum(1 for r in results if r[2] == 'FAIL')
skipped = sum(1 for r in results if r[2] == 'SKIP')

print(f'\n=== {nb_path} SUMMARY: {passed} PASS, {failed} FAIL, {skipped} SKIP ===')
for r in results:
    status = 'PASS' if r[2] == 'PASS' else ('FAIL' if r[2] == 'FAIL' else 'SKIP')
    print(f'  [{status}] cell[{r[1]}] #{r[0]}: {r[3][:120]}')

if failed > 0:
    sys.exit(1)
