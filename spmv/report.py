"""Summarise CSR SpMV abstraction results: ratio = t_variant / t_f77.

Reads every results/*.txt written by run.sh; the column label comes from the
'# compiler:' and '# flags:' header lines of each file.
"""
import glob
import math
import re

gm = lambda v: math.exp(sum(map(math.log, v)) / len(v))


def load(path):
    compiler, flags, rows = '?', '', []
    for line in open(path):
        if line.startswith('# compiler:'):
            compiler = line.split(':', 1)[1].strip()
        elif line.startswith('# flags:'):
            flags = line.split(':', 1)[1].strip()
        elif not line.startswith('#') and line.strip():
            m, v, n, nnz, tf, tv = line.split()
            rows.append((m, v, int(n), int(nnz), float(tf), float(tv)))
    opt = re.search(r'-O\d', flags)
    label = f"{compiler[:7]} {opt.group(0) if opt else ''}"
    if 'native' in flags:
        label += 'n'
    return label, rows


runs = sorted((load(p) for p in glob.glob('results/*.txt')
               if not p.endswith('summary.txt')), key=lambda r: r[0])
if not runs:
    raise SystemExit('no results/*.txt found; run ./run.sh first')
labels = [r[0] for r in runs]
data = [r[1] for r in runs]
mats = list(dict.fromkeys(r[0] for r in data[0]))
vars_ = list(dict.fromkeys(r[1] for r in data[0]))


def row(label, vals, fmt='{:13.2f}'):
    print(f'{label:24s}' + ''.join(fmt.format(v) for v in vals))


def header():
    print(f"{'':24s}" + ''.join(f'{s:>13s}' for s in labels))


print('CSR SpMV: ratio t_variant / t_f77, 1.00 = no abstraction penalty.')
print("'n' after the optimisation level = native CPU tuning.")
print('ctrl is the F77 loop copied verbatim: its spread is the noise floor.\n')
print('Matrices:')
for r in data[0]:
    if r[1] == vars_[0]:
        print(f'  {r[0]:12s} n={r[2]:<8d} nnz={r[3]:<9d} ({r[3]/r[2]:.1f} per row)')

print('\nGeomean over all variants and matrices\n'); header()
row('all', [gm([r[5]/r[4] for r in d]) for d in data])
for m in mats:
    row(m, [gm([r[5]/r[4] for r in d if r[0] == m]) for d in data])

print('\nPer variant, geomean over matrices\n'); header()
for v in vars_:
    row(v, [gm([r[5]/r[4] for r in d if r[1] == v]) for d in data])

for m in mats:
    print(f'\nPer variant, {m}\n'); header()
    for v in vars_:
        row(v, [next(r[5]/r[4] for r in d if r[0] == m and r[1] == v) for d in data])

print('\nF77 baseline, ns per nonzero\n'); header()
for m in mats:
    row(m, [next(r[4]/r[3] for r in d if r[0] == m) for d in data], '{:13.3f}')
