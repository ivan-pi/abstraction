import glob, math, collections
rows = collections.OrderedDict()
order = ["gfortran_O2","gfortran_O3marchnative","gfortran_O3marchnativeffastmath",
         "flang19_O2","flang19_O3marchnative","flang19_O3marchnativeffastmath"]
data = {}
for cid in order:
    lines = open(f"results/{cid}.txt").read().splitlines()
    tag = lines[0].split(": ")[1] + " " + lines[1].split(": ")[1]
    d = []
    for l in lines:
        if l.startswith("#"): continue
        k, n, tl, ta = l.split()
        d.append((k, int(n), float(tl), float(ta)))
    data[tag] = d
gm = lambda v: math.exp(sum(map(math.log, v))/len(v))
sizes = sorted({r[1] for d in data.values() for r in d})
kernels = list(dict.fromkeys(r[0] for d in data.values() for r in d))
print("Geometric mean of t_array / t_loop   (>1 = array syntax slower)\n")
print(f"{'config':40s}" + "".join(f"{'n='+str(s):>11s}" for s in sizes) + f"{'all':>9s}")
for tag, d in data.items():
    per = [gm([r[3]/r[2] for r in d if r[1]==s]) for s in sizes]
    print(f"{tag:40s}" + "".join(f"{p:11.3f}" for p in per) + f"{gm([r[3]/r[2] for r in d]):9.3f}")
print("\nPer-kernel ratio t_array/t_loop at n=1000 (L1-resident)\n")
tags = list(data)
short = ["gf -O2","gf -O3 nat","gf fast","fl -O2","fl -O3 nat","fl fast"]
print(f"{'kernel':9s}" + "".join(f"{s:>11s}" for s in short))
for k in kernels:
    vals = [next(r[3]/r[2] for r in data[t] if r[0]==k and r[1]==1000) for t in tags]
    print(f"{k:9s}" + "".join(f"{v:11.2f}" for v in vals))
print("\nAbsolute loop-version time at n=1000, ns/call\n")
print(f"{'kernel':9s}" + "".join(f"{s:>11s}" for s in short))
for k in kernels:
    vals = [next(r[2] for r in data[t] if r[0]==k and r[1]==1000) for t in tags]
    print(f"{k:9s}" + "".join(f"{v:11.0f}" for v in vals))
for s in sizes[1:]:
    print(f"\nPer-kernel ratio at n={s}\n")
    print(f"{'kernel':9s}" + "".join(f"{x:>11s}" for x in short))
    for k in kernels:
        vals = [next(r[3]/r[2] for r in data[t] if r[0]==k and r[1]==s) for t in tags]
        print(f"{k:9s}" + "".join(f"{v:11.2f}" for v in vals))
