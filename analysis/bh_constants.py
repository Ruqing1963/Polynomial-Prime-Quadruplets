#!/usr/bin/env python3
"""
bh_constants.py -- Bateman-Horn singular series C_k for the k-tuple Q(n), Q(n+1), ..., Q(n+k-1),
Q(n) = n^q - (n-1)^q, and the predicted k-tuple counts up to Y.

  C_k = prod_p (1 - w_k(p)/p) / (1 - 1/p)^k,   w_k(p) = #{n mod p : p | Q(n)...Q(n+k-1)}.

Q has no roots mod p unless p == 1 (mod q), and then exactly q-1 roots, n = 1/(1-u) with u a
nontrivial q-th root of unity.  So only primes p == 1 (mod q) carry w_k > 0; every other prime
contributes (1 - 1/p)^-k.  The product converges only conditionally (the mean of w_1 is 1);
we report partial products at several cutoffs M to show the convergence.

Usage:  python bh_constants.py --qs 23 29 37 41 43 47 53 59 61 71 --pmax 2e8 --out bh_constants.json
"""
import argparse
import json
import math
from multiprocessing import Pool

import numpy as np


def primes_upto(m):
    s = np.ones(m + 1, dtype=bool)
    s[:2] = False
    s[4::2] = False
    for i in range(3, math.isqrt(m) + 1, 2):
        if s[i]:
            s[i * i :: 2 * i] = False
    return np.nonzero(s)[0]


def roots(p, q):
    """The q-1 residues n mod p with p | Q(n), for a prime p == 1 (mod q)."""
    e = (p - 1) // q
    g = 2
    while (w := pow(g, e, p)) == 1:
        g += 1
    out, u = [], 1
    for _ in range(q - 1):
        u = u * w % p
        out.append(pow(1 - u, -1, p))
    return out


def run(args):
    q, ps, cuts, kmax = args
    logsum = np.zeros(kmax + 1)  # logsum[k] = log of partial product of C_k
    res = {k: {} for k in range(1, kmax + 1)}
    ci = 0
    for p in ps:
        p = int(p)
        while ci < len(cuts) and p > cuts[ci]:
            for k in range(1, kmax + 1):
                res[k][cuts[ci]] = math.exp(logsum[k])
            ci += 1
        base = -math.log1p(-1.0 / p)
        if p % q == 1:
            r = roots(p, q)
            s = set()
            for k in range(1, kmax + 1):
                s.update((x - (k - 1)) % p for x in r)
                logsum[k] += math.log1p(-len(s) / p) + k * base
        else:
            for k in range(1, kmax + 1):
                logsum[k] += k * base
    for c in cuts[ci:]:
        for k in range(1, kmax + 1):
            res[k][c] = math.exp(logsum[k])
    return q, res


def expected(Ck, q, k, Y, y0=2.0):
    """C_k * int_{y0}^{Y} dt / prod_i ln Q(t+i),  ln Q(t) ~ ln q + (q-1) ln t."""
    x = np.linspace(math.log(y0), math.log(Y), 200001)
    t = np.exp(x)
    f = t / (math.log(q) + (q - 1) * x) ** k
    return Ck * float(np.sum((f[1:] + f[:-1]) / 2 * np.diff(x)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qs", type=int, nargs="+", default=[23, 29, 37, 41, 43, 47, 53, 59, 61, 71])
    ap.add_argument("--pmax", type=float, default=2e8)
    ap.add_argument("--kmax", type=int, default=4)
    ap.add_argument("--out", default="bh_constants.json")
    a = ap.parse_args()
    M = int(a.pmax)
    cuts = sorted({int(M / 100), int(M / 10), M})
    ps = primes_upto(M)
    with Pool(min(len(a.qs), 8)) as pool:
        results = dict(pool.map(run, [(q, ps, cuts, a.kmax) for q in a.qs]))
    out = {}
    for q in a.qs:
        out[q] = {}
        for k in range(1, a.kmax + 1):
            C = results[q][k][M]
            out[q][k] = {
                "C_partial": {str(c): results[q][k][c] for c in cuts},
                "E(1e8)": expected(C, q, k, 1e8),
                "E(2e9)": expected(C, q, k, 2e9),
            }
        row = "  ".join(
            f"C{k}={out[q][k]['C_partial'][str(M)]:.4g} E4(1e8)={out[q][k]['E(1e8)']:.3g}" if k == 4 else
            f"C{k}={out[q][k]['C_partial'][str(M)]:.4g}" for k in range(1, a.kmax + 1))
        print(f"q={q:3d}  {row}  E4(2e9)={out[q][a.kmax]['E(2e9)']:.3g}", flush=True)
    with open(a.out, "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
