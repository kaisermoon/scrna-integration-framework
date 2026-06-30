#!/usr/bin/env python3
"""Execute 07_deg.ipynb cells in order, reporting pass/fail for each."""

import nbformat, sys, os, traceback, time
from pathlib import Path

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
os.chdir(_PROJECT_ROOT)
nb_path = "notebooks/07_downstream/07_deg.ipynb"
os.makedirs('results/figures', exist_ok=True)
os.makedirs('results/tables', exist_ok=True)

nb = nbformat.read(nb_path, as_version=4)
code_cells = [(i, c) for i, c in enumerate(nb.cells) if c.cell_type == 'code']

print(f'=== {nb_path} ({len(code_cells)} code cells) ===')

results = []
for idx, (cell_idx, cell) in enumerate(code_cells):
    cell_num = idx + 1
    src = cell.source
    if not src.strip():
        results.append((cell_num, cell_idx, 'SKIP', 'empty'))
        continue

    print(f'[{cell_num}/{len(code_cells)}] cell[{cell_idx}]...', end=' ', flush=True)
    t0 = time.time()
    try:
        exec(src, globals())
        elapsed = time.time() - t0
        print(f'OK ({elapsed:.1f}s)')
        results.append((cell_num, cell_idx, 'PASS', f'{elapsed:.1f}s'))
    except Exception as e:
        elapsed = time.time() - t0
        err = traceback.format_exc()
        print(f'FAIL ({elapsed:.1f}s)')
        lines = err.strip().split('\n')
        for l in lines[-8:]:
            print(f'  {l}')
        results.append((cell_num, cell_idx, 'FAIL', str(e)[:200]))

# Summary
passed = sum(1 for r in results if r[2] == 'PASS')
failed = sum(1 for r in results if r[2] == 'FAIL')
skipped = sum(1 for r in results if r[2] == 'SKIP')

print(f'\n=== {nb_path} SUMMARY: {passed} PASS, {failed} FAIL, {skipped} SKIP ===')
for r in results:
    if r[2] == 'PASS':
        status = '✓'
    elif r[2] == 'FAIL':
        status = '✗'
    else:
        status = '-'
    print(f'  {status} cell[{r[1]}] #{r[0]}: {r[3][:120]}')

# Check for output files
output_h5ad = "results/07_deg_v1.h5ad"
if os.path.exists(output_h5ad):
    print(f'\nOutput file produced: {output_h5ad} ({os.path.getsize(output_h5ad):,} bytes)')
else:
    print(f'\nOutput file NOT produced: {output_h5ad}')

import glob
figs = glob.glob('results/figures/07_deg_*.png')
if figs:
    print(f'Figure files ({len(figs)}):')
    for f in sorted(figs):
        print(f'  {f} ({os.path.getsize(f):,} bytes)')
else:
    print('No figure files produced.')

tables = glob.glob('results/tables/07_deg_*.csv') + glob.glob('results/tables/07_deg_*.tsv')
if tables:
    print(f'Table files ({len(tables)}):')
    for f in sorted(tables):
        print(f'  {f} ({os.path.getsize(f):,} bytes)')
else:
    print('No table files produced.')

if failed > 0:
    sys.exit(1)
