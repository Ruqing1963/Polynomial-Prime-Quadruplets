#!/usr/bin/env python3
"""
patch_prime_gaps.py -- re-test every gap of more than GAP consecutive n with no listed Q47 prime
(from count_primes_2e9.py's bitmap) and add any primes the data files missed.

A gap > 2000 near n = 1e9 has probability about exp(-2000/105) ~ 5e-9 per prime, so over 1.8e7 primes
fewer than one is expected by chance; each one is tested directly (strong base-2 PRP, then BPSW).
Output: patch_prime_gaps.json, with corrected counts of primes, doublets and triplets for n <= 2e9.

Usage:  python patch_prime_gaps.py <Data_n_1_to_300M.zip or directory> <Data_n=300M_plus directory>
"""
import json
import sys

import gmpy2
import numpy as np

import count_primes_2e9 as cp

GAP = 2000


def main():
    sys.argv = sys.argv[:3]
    zpath, ddir = sys.argv[1], sys.argv[2]
    bits = np.zeros(cp.Y + 2, dtype=bool)
    import glob, os, zipfile
    if os.path.isdir(zpath):
        blobs = [open(f, "rb").read() for f in glob.glob(os.path.join(zpath, "*.txt"))]
    else:
        with zipfile.ZipFile(zpath) as z:
            blobs = [z.read(e) for e in z.infolist() if e.filename.endswith(".txt")]
    for data in blobs:
        a = cp.ints(data)
        bits[a[(a >= 1) & (a <= cp.Y)]] = True
    del blobs
    for f in glob.glob(os.path.join(ddir, "Prime_Q47_*M_*M.txt")):
        a = cp.ints(open(f, "rb").read())
        bits[a[(a >= 1) & (a <= cp.Y)]] = True
    listed = int(bits.sum())
    n = np.flatnonzero(bits)
    gaps = [(int(n[i]), int(n[i + 1])) for i in np.flatnonzero(np.diff(n) > GAP)]
    added = {}
    for lo, hi in gaps:
        new = [m for m in range(lo + 1, hi)
               if gmpy2.is_strong_prp(gmpy2.mpz(m) ** 47 - gmpy2.mpz(m - 1) ** 47, 2)
               and gmpy2.is_bpsw_prp(gmpy2.mpz(m) ** 47 - gmpy2.mpz(m - 1) ** 47)]
        added[f"{lo}-{hi}"] = new
        bits[new] = True
        print(lo, hi, len(new), flush=True)
    b = bits[1:cp.Y + 1]
    out = {"listed": listed, "gaps_over_%d" % GAP: gaps, "added": added,
           "primes": int(b.sum()), "doublets": int(np.sum(b[:-1] & b[1:])),
           "triplets": int(np.sum(b[:-2] & b[1:-1] & b[2:])),
           "quadruplets": int(np.sum(b[:-3] & b[1:-2] & b[2:-1] & b[3:]))}
    print(json.dumps({k: v for k, v in out.items() if k != "added"}, indent=1))
    json.dump(out, open("patch_prime_gaps.json", "w"), indent=1)


if __name__ == "__main__":
    main()
