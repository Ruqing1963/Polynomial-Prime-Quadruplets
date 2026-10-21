#!/usr/bin/env python3
"""
window_1e44.py -- summarize the author's Q47 prime lists in two windows near 10^44 and compare with
Bateman-Horn: expected count = C1 * sum_n 1/ln Q47(n) over the window (ln Q47(n) ~ ln 47 + 46 ln n).
Also counts doublets/triplets inside each list and checks a random sample with BPSW.
Output: window_1e44.json.

Usage:  python window_1e44.py <file> [<file> ...]
"""
import json
import math
import random
import re
import sys

import gmpy2

C1 = json.load(open("bh_constants.json"))["47"]["1"]["C_partial"]["200000000"]
C2 = json.load(open("bh_constants.json"))["47"]["2"]["C_partial"]["200000000"]
out = {}
for path in sys.argv[1:]:
    ns = sorted({int(m) for m in re.findall(r"^\s*(\d{20,})", open(path).read(), re.M)})
    lo, hi = ns[0], ns[-1]
    # window bounds: take the file name's range if it states one, else the min/max of the data
    L = hi - lo + 1
    lnQ = math.log(47) + 46 * math.log(lo + L / 2)
    S = set(ns)
    dbl = sum(1 for n in ns if n + 1 in S)
    tri = sum(1 for n in ns if n + 1 in S and n + 2 in S)
    sample = random.Random(1).sample(ns, min(40, len(ns)))
    ok = sum(gmpy2.is_bpsw_prp(gmpy2.mpz(n) ** 47 - gmpy2.mpz(n - 1) ** 47) for n in sample)
    out[path.split("\\")[-1]] = {
        "count": len(ns), "min_n": str(lo), "max_n": str(hi), "span": L,
        "expected_BH_over_span": C1 * L / lnQ, "doublets": dbl,
        "expected_doublets": C2 * L / lnQ ** 2, "triplets": tri,
        "bpsw_sample": f"{ok}/{len(sample)}"}
    print(json.dumps(out[path.split("\\")[-1]]), flush=True)
json.dump(out, open("window_1e44.json", "w"), indent=1)
