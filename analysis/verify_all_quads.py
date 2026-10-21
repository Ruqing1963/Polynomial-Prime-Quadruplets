#!/usr/bin/env python3
"""
verify_all_quads.py -- BPSW check of every Q47 quadruplet start n < 8e11 (zip census + the 15 below 2e9).
Output: verify_all_quads.json.
"""
import json
from multiprocessing import Pool

import gmpy2

V1 = [23159557, 117309848, 136584738, 218787064, 411784485, 423600750, 523331634, 640399031,
      987980498, 1163461515, 1370439187, 1643105964, 1691581855, 1975860550, 1996430175]


def Q(n):
    n = gmpy2.mpz(n)
    return n ** 47 - (n - 1) ** 47


def check(n):
    return n, all(gmpy2.is_bpsw_prp(Q(n + j)) for j in range(4))


if __name__ == "__main__":
    S = sorted({n for n in json.load(open("quad_starts_from_zip.json")) if n < 8e11} | set(V1))
    with Pool(8) as pool:
        res = dict(pool.map(check, S, chunksize=16))
    bad = [n for n, ok in res.items() if not ok]
    good = sorted(n for n, ok in res.items() if ok)
    out = {"tested": len(S), "bpsw_ok": len(good), "failed": bad,
           "count_lt_2e9": sum(n < 2e9 for n in good),
           "count_lt_2e11": sum(n < 2e11 for n in good),
           "count_2.3e7_8e11": sum(2.3e7 <= n < 8e11 for n in good),
           "count_lt_8e11": len(good)}
    print(json.dumps(out, indent=1))
    json.dump(out, open("verify_all_quads.json", "w"), indent=1)
