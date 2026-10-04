"""Summarise BLAS-2 abstraction results: ratio = t_variant / t_f77."""
import glob, math, collections

FACILITY = {
    'arrsec': 'array sections', 'arrsec_tmp': 'array sections', 'arrsec_sum': 'array sections',
    'ashape': 'assumed shape',
    'matmul': 'intrinsics', 'sum': 'intrinsics', 'dotprod': 'intrinsics',
    'matmul_vm': 'intrinsics', 'matmul_tr': 'intrinsics', 'spread': 'intrinsics',
    'dc_ji': 'do concurrent', 'dc_ij': 'do concurrent',
    'dtype': 'derived types', 'tbp': 'derived types', 'pdt': 'derived types',
    'accessor': 'derived types',
}
ORDER = ['gfortran_O2', 'gfortran_O3marchnative', 'flang19_O2', 'flang19_O3marchnative']
SHORT = ['gf -O2', 'gf -O3 nat', 'fl -O2', 'fl -O3 nat']
gm = lambda v: math.exp(sum(map(math.log, v)) / len(v))

data = {}
nonbitwise = {}   # (config) -> set of (kernel, variant) whose result differs from F77
for cid in ORDER:
    rows = []
    nonbitwise[cid] = set()
    for line in open(f'results/{cid}.txt'):
        if line.startswith('# note: not bitwise:'):
            f = line.split(':')[2].split()
            nonbitwise[cid].add((f[0], f[1]))
        if line.startswith('#'):
            continue
        k, v, n, tf, tv = line.split()
        rows.append((k, v, int(n), float(tf), float(tv)))
    data[cid] = rows
sizes = sorted({r[2] for r in data[ORDER[0]]})
pairs = list(dict.fromkeys((r[0], r[1]) for r in data[ORDER[0]]))

def ratio(cid, k, v, n):
    for r in data[cid]:
        if r[:3] == (k, v, n):
            return r[4] / r[3]
    return None

def table(title, keys, fn):
    print(f'\n{title}\n')
    print(f"{'':22s}" + ''.join(f'{s:>12s}' for s in SHORT))
    for key in keys:
        vals = [fn(cid, key) for cid in ORDER]
        print(f'{key:22s}' + ''.join(f'{v:12.2f}' if v else f"{'n/a':>12s}" for v in vals))

print('Ratio t_variant / t_f77, geometric means.  1.00 = no abstraction penalty.')
print('Geomeans exclude variants whose result is not bitwise identical to F77')
print('(they did different arithmetic):')
for cid, sh in zip(ORDER, SHORT):
    print(f'  {sh:12s}', ', '.join(sorted(f'{k}/{v}' for k, v in nonbitwise[cid])) or 'none')

ok = lambda cid, r: (r[0], r[1]) not in nonbitwise[cid]
table('Overall, per size', [f'n={n}' for n in sizes] + ['all sizes'],
      lambda cid, key: gm([r[4]/r[3] for r in data[cid] if ok(cid, r)
                           and (key == 'all sizes' or r[2] == int(key[2:]))]))

facs = list(dict.fromkeys(FACILITY[v] for _, v in pairs))
table('Per language facility, all kernels and sizes', facs,
      lambda cid, f: gm([r[4]/r[3] for r in data[cid] if ok(cid, r) and FACILITY[r[1]] == f])
      if any(ok(cid, r) and FACILITY[r[1]] == f for r in data[cid]) else None)

table('Per language facility, excluding ger/spread', facs,
      lambda cid, f: gm([r[4]/r[3] for r in data[cid] if ok(cid, r) and FACILITY[r[1]] == f
                         and r[1] != 'spread'])
      if any(ok(cid, r) and FACILITY[r[1]] == f for r in data[cid]) else None)

table('Per variant, geometric mean over all sizes  (* = not bitwise for some compiler)',
      [f'{k}/{v}' + ('*' if any((k, v) in nonbitwise[c] for c in ORDER) else '') for k, v in pairs],
      lambda cid, key: gm([r[4]/r[3] for r in data[cid] if f'{r[0]}/{r[1]}' == key.rstrip('*')])
      if any(f'{r[0]}/{r[1]}' == key.rstrip('*') for r in data[cid]) else None)

for n in sizes:
    table(f'Per variant, n={n}', [f'{k}/{v}' for k, v in pairs],
          lambda cid, key: ratio(cid, *key.split('/'), n))

print('\nF77 baseline time, ns/call\n')
print(f"{'':22s}" + ''.join(f'{s:>12s}' for s in SHORT))
for k in dict.fromkeys(p[0] for p in pairs):
    for n in sizes:
        vals = [next(r[3] for r in data[cid] if r[0] == k and r[2] == n) for cid in ORDER]
        print(f'{k+" n="+str(n):22s}' + ''.join(f'{v:12.0f}' for v in vals))
