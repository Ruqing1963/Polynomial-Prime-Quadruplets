#!/usr/bin/env python3
"""
control_search.py -- count primes and k-tuples (k = 2, 3, 4) of Q(n) = n^q - (n-1)^q for 2 <= n <= Y,
for several exponents q, and compare with the Bateman-Horn predictions from bh_constants.json.

Each n is first sieved with the primes p == 1 (mod q), p <= B (all other primes never divide Q(n));
survivors get one strong probable-prime test to base 2 (gmpy2).  Q(n) > B for every n >= 2, so a
sieve hit always means "composite".

Usage:  python control_search.py --qs 23 29 41 43 47 53 59 71 --Y 1e7 --workers 16
"""
import argparse
import json
import math
import time
from multiprocessing import Pool

import gmpy2
import numpy as np

from bh_constants import expected, primes_upto, roots

CHUNK = 250_000


def sieve_chunk(args):
    q, lo, hi, sp = args  # n in [lo, hi)
    alive = np.ones(hi - lo, dtype=bool)
    for p in sp:
        for r in roots(p, q):
            alive[(r - lo) % p :: p] = False
    out = np.zeros(hi - lo, dtype=bool)
    for i in np.nonzero(alive)[0].tolist():
        n = gmpy2.mpz(lo + i)
        out[i] = gmpy2.is_strong_prp(n ** q - (n - 1) ** q, 2)
    return lo, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qs", type=int, nargs="+", default=[23, 29, 41, 43, 47, 53, 59, 71])
    ap.add_argument("--Y", type=float, default=1e7)
    ap.add_argument("--B", type=float, default=1e5)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--consts", default="bh_constants.json")
    ap.add_argument("--out", default="control_search.json")
    a = ap.parse_args()
    Y, B = int(a.Y), int(a.B)
    consts = json.load(open(a.consts))
    P = primes_upto(B)
    results = {}
    with Pool(a.workers) as pool:
        for q in a.qs:
            t0 = time.time()
            sp = [int(p) for p in P if p % q == 1]
            jobs = [(q, lo, min(lo + CHUNK, Y + 1), sp) for lo in range(2, Y + 1, CHUNK)]
            isp = np.zeros(Y + 4, dtype=bool)
            for lo, out in pool.imap_unordered(sieve_chunk, jobs):
                isp[lo : lo + len(out)] = out
            row = {"Y": Y, "time_s": round(time.time() - t0)}
            run = isp[2 : Y + 1].astype(np.int8)
            for k in range(1, 5):
                obs = int(np.sum(np.all([run[i : len(run) - k + 1 + i] for i in range(k)], axis=0)))
                C = consts[str(q)][str(k)]["C_partial"]["200000000"]
                E = expected(C, q, k, Y)
                row[k] = {"observed": obs, "expected": E, "O/E": obs / E if E else None,
                          "z": (obs - E) / math.sqrt(E) if E else None}
            results[q] = row
            print(q, json.dumps(row), flush=True)
            json.dump(results, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
