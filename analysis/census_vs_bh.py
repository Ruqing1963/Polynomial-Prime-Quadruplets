import json, math
from bh_constants import expected
C = json.load(open("bh_constants_q47_k6.json"))["47"]
Ck = {k: C[str(k)]["C_partial"]["200000000"] for k in range(1, 7)}
print({k: round(v, 1) for k, v in Ck.items()})
S = sorted(set(json.load(open("quad_starts_from_zip.json"))))
paper15 = [23159557,117309848,136584738,218787064,411784485,423600750,523331634,640399031,987980498,1163461515,1370439187,1643105964,1691581855,1975860550,1996430175]
low = sorted(set(S) | set(paper15))
cen = [n for n in low if n < 8e11]
print("census <8e11:", len(cen))
E = lambda k, a, b: expected(Ck[k], 47, k, b) - expected(Ck[k], 47, k, a)
edges = [2, 2e9, 1e10, 3e10, 1e11, 2e11, 4e11, 6e11, 8e11]
tot_o = tot_e = 0
for a, b in zip(edges, edges[1:]):
    o = sum(a <= n < b for n in cen); e = E(4, a, b); tot_o += o; tot_e += e
    print(f"[{a:.0e},{b:.0e}) obs {o:5d}  BH {e:8.1f}  O/E {o/e:.3f}  z {(o-e)/math.sqrt(e):+.2f}")
print(f"total obs {tot_o} BH {tot_e:.1f} O/E {tot_o/tot_e:.4f} z {(tot_o-tot_e)/math.sqrt(tot_e):+.2f}")
st = set(cen)
q5 = [n for n in cen if n + 1 in st]; q6 = [n for n in cen if n + 1 in st and n + 2 in st]
print("quintuplets <8e11:", len(q5), q5, " BH E5:", round(E(5, 2, 8e11), 2), " sextuplets:", len(q6), "BH E6:", round(E(6, 2, 8e11), 3))
# deep windows (end = last hit, so E is a lower bound on the window's expectation)
for name, a, last, cnt in [("1e20 window", 1e20, 100000000030413600025, 14),
                           ("1000-digit window", 4797791281515011781318, 4797791281552013327611, 11)]:
    L = last - a; lnQ = math.log(47) + 46 * math.log(a)
    e = Ck[4] * L / lnQ ** 4
    print(f"{name}: length >= {L:.3e}, digits {lnQ/math.log(10):.0f}, obs {cnt}, BH {e:.2f} (to last hit)")
