import glob, re, json, os, math
from collections import Counter
root = r"zip_4tuplets"
starts = set(); per_file = {}
for f in glob.glob(root + r"/**/*", recursive=True):
    if not os.path.isfile(f) or not f.endswith((".log", ".txt")): continue
    s = open(f, encoding="utf-8", errors="ignore").read()
    m = [int(x) for x in re.findall(r"(?:Sequence|QUADRUPLET):\s*(\d+)", s)]
    if m:
        per_file[os.path.relpath(f, root)] = (len(m), min(m), max(m))
        starts.update(m)
for k, v in sorted(per_file.items(), key=lambda kv: kv[1][1]): print(v, k)
S = sorted(starts)
print("distinct quadruplet starts:", len(S), "min", S[0], "max", S[-1])
edges = [0, 2e9, 4e9, 8e9, 1e10, 1.5e10, 2e10, 3e10, 1e11, 1.5e11, 2e11, 3e11, 4e11, 5e11, 6e11, 7e11, 8e11]
c = Counter()
for n in S:
    for a, b in zip(edges, edges[1:]):
        if a <= n < b: c[(a, b)] += 1
for a, b in zip(edges, edges[1:]): print(f"[{a:.3g},{b:.3g}) {c[(a,b)]}")
json.dump(S, open("quad_starts_from_zip.json", "w"))
