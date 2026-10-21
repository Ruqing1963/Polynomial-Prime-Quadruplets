#!/usr/bin/env python3
"""
rescan_window.py -- independent complete rescan of a window lo <= n < hi for primes Q47(n), and comparison
with an existing prime list.

Method differs from the original search (Miller-Rabin on every n): every n is first sieved with the primes
p == 1 (mod 47), p <= B (only these can divide Q47(n), on 46 classes each; Q47(n) > B so a hit proves
compositeness); survivors get a Fermat base-3 PRP test with PFGW (one process per core, ABC input
files), and every PFGW hit is confirmed with BPSW (gmpy2).

Usage:  python rescan_window.py --lo <int> --hi <int> --list <prime list file> --tag w1
                                --pfgw <path to pfgw64.exe> --work <empty work dir> [--B 1e8]
Output: rescan_<tag>.json
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import time
from multiprocessing import Pool

import gmpy2
import numpy as np

from bh_constants import primes_upto, roots


def bpsw(m):
    n = gmpy2.mpz(m)
    return m, gmpy2.is_bpsw_prp(n ** 47 - (n - 1) ** 47)


def run_pfgw(exe, work, ns, k):
    """Split ns over k PFGW processes; return the set of n whose Q47(n) PFGW reports as PRP."""
    procs = []
    for j in range(k):
        d = os.path.join(work, f"w{j:02d}")
        os.makedirs(d, exist_ok=True)
        shutil.copy(exe, d)
        with open(os.path.join(d, "in.txt"), "w") as f:
            f.write("ABC ($a)^47-($a-1)^47\n")
            f.writelines(f"{n}\n" for n in ns[j::k])
        procs.append(subprocess.Popen([os.path.join(d, "pfgw64.exe"), "-f0", "-k", "in.txt"],
                                      cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL))
    for p in procs:
        p.wait()
    hits = set()
    for j in range(k):
        log = os.path.join(work, f"w{j:02d}", "pfgw.log")
        if os.path.exists(log):
            hits |= {int(m) for m in re.findall(r"^\((\d+)\)\^47", open(log).read(), re.M)}
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lo", type=int, required=True)
    ap.add_argument("--hi", type=int, required=True)
    ap.add_argument("--list", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--B", type=float, default=1e8)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--pfgw", required=True)
    ap.add_argument("--work", required=True)
    a = ap.parse_args()
    t0, L = time.time(), a.hi - a.lo
    alive = np.ones(L, dtype=bool)
    for p in primes_upto(int(a.B)):
        p = int(p)
        if p % 47 != 1:
            continue
        o = (np.array(roots(p, 47), dtype=np.int64) - a.lo % p) % p
        if p < L:
            for x in o.tolist():
                alive[x::p] = False
        else:
            alive[o[o < L]] = False
    surv = (np.flatnonzero(alive) + 0).tolist()
    print(f"{a.tag}: sieve done, {len(surv)} survivors of {L} [{time.time() - t0:.0f}s]", flush=True)
    tests = len(surv)
    hits = run_pfgw(a.pfgw, a.work, [a.lo + i for i in surv], a.workers)
    print(f"{a.tag}: PFGW done, {len(hits)} PRPs [{time.time() - t0:.0f}s]", flush=True)
    with Pool(a.workers) as pool:
        res = pool.map(bpsw, sorted(hits), chunksize=8)
    found = sorted(m for m, ok in res if ok)
    print(f"{a.tag}: BPSW confirmed {len(found)} [{time.time() - t0:.0f}s]", flush=True)
    old = {int(m) for m in re.findall(r"^\s*(\d{20,})", open(a.list, encoding="utf-8").read(), re.M)}
    old = {n for n in old if a.lo <= n < a.hi}
    new = set(found)
    out = {"lo": str(a.lo), "hi": str(a.hi), "B": a.B, "primes": len(found), "prp_tests": tests,
           "pfgw_prp_not_bpsw": sorted(str(m) for m, ok in res if not ok),
           "listed": len(old), "missing_from_list": sorted(str(n) for n in new - old),
           "not_confirmed": sorted(str(n) for n in old - new),
           "doublets": sum(1 for n in found if n + 1 in new),
           "triplets": sum(1 for n in found if n + 1 in new and n + 2 in new),
           "time_s": round(time.time() - t0)}
    print(json.dumps(out, indent=1), flush=True)
    json.dump(out, open(f"rescan_{a.tag}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
