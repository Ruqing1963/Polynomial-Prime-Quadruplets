#!/usr/bin/env python3
"""
quad_census.py -- complete census of prime quadruplets Q(n), Q(n+1), Q(n+2), Q(n+3), 2 <= n <= Y,
for Q(n) = n^q - (n-1)^q.

Method: a segmented sieve marks every n for which Q(n) has a prime factor p <= B (only p == 1 mod q
can divide Q(n), on q-1 residue classes each).  An n survives as a quadruplet candidate only if
n, n+1, n+2, n+3 all survive.  Candidates get strong base-2 probable-prime tests on Q(n), ..., Q(n+3)
with early exit; every quadruplet found is re-checked with BPSW.  Since Q(n) >= 2^q - 1 > B for
n >= 2, a sieve hit always proves compositeness, so the census is exact up to PRP reliability.

Usage:  python quad_census.py --q 47 --Y 2e9 --workers 14
"""
import argparse
import json
import math
import time
from multiprocessing import Pool

import gmpy2
import numpy as np

from bh_constants import expected, primes_upto, roots

_W = {}


def _init(q, table):
    _W["q"], _W["table"] = q, table


def segment(bounds):
    lo, hi = bounds  # candidate starts n in [lo, hi); sieve covers [lo, hi + 3)
    q = _W["q"]
    L = hi + 3 - lo
    alive = np.ones(L, dtype=bool)
    for p, rs in _W["table"]:
        for o in ((rs - lo) % p).tolist():
            alive[o::p] = False
    quad = alive[:-3] & alive[1:-2] & alive[2:-1] & alive[3:]
    cand = np.nonzero(quad)[0]
    found, tests, cache = [], 0, {}
    for i in cand.tolist():
        n0 = lo + i
        ok = True
        for j in range(4):
            m = n0 + j
            if m not in cache:
                n = gmpy2.mpz(m)
                tests += 1
                cache[m] = gmpy2.is_strong_prp(n ** q - (n - 1) ** q, 2)
            if not cache[m]:
                ok = False
                break
        if ok:
            found.append(n0)
    return lo, hi, len(cand), tests, found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--Y", type=float, default=2e9)
    ap.add_argument("--B", type=float, default=3e7)
    ap.add_argument("--seg", type=float, default=1e7)
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--consts", default="bh_constants.json")
    a = ap.parse_args()
    q, Y, B, S = a.q, int(a.Y), int(a.B), int(a.seg)
    table = [(int(p), np.array(roots(int(p), q), dtype=np.int64))
             for p in primes_upto(B) if p % q == 1]
    jobs = [(lo, min(lo + S, Y + 1)) for lo in range(2, Y + 1, S)]
    t0 = time.time()
    found, ncand, ntest = [], 0, 0
    with Pool(a.workers, _init, (q, table)) as pool:
        for i, (lo, hi, c, t, f) in enumerate(pool.imap_unordered(segment, jobs), 1):
            found += f
            ncand += c
            ntest += t
            if i % 10 == 0 or i == len(jobs):
                print(f"q={q} {i}/{len(jobs)} segments, {len(found)} quadruplets, "
                      f"{ncand} candidates, {ntest} PRP tests [{time.time() - t0:.0f}s]", flush=True)
    found.sort()
    verified = [n for n in found
                if all(gmpy2.is_bpsw_prp(gmpy2.mpz(n + j) ** q - gmpy2.mpz(n + j - 1) ** q)
                       for j in range(4))]
    consts = json.load(open(a.consts))
    E = expected(consts[str(q)]["4"]["C_partial"]["200000000"], q, 4, Y)
    out = {"q": q, "Y": Y, "B": B, "quadruplets": verified, "count": len(verified),
           "expected_BH": E, "candidates": ncand, "prp_tests": ntest,
           "time_s": round(time.time() - t0)}
    print(json.dumps(out), flush=True)
    json.dump(out, open(f"quad_census_q{q}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
