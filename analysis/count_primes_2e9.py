#!/usr/bin/env python3
"""
count_primes_2e9.py -- count the Q47 prime list (n with Q47(n) prime, 1 <= n <= 2e9) from the author's
data files: Data_n_1_to_300M.zip and Data_n=300M_plus/Prime_Q47_*M_*M.txt.

Each line starting with an integer n (or "n = <n>") is one prime.  The script builds a bitmap over n, reports counts per
file, duplicates, coverage gaps (runs of > GAP consecutive n with no prime, which have probability
~ exp(-GAP/104) near n = 1e9 and so indicate missing data), and the numbers of primes, doublets and
triplets for n <= 2e9.  Output: count_primes_2e9.json.

Usage:  python count_primes_2e9.py <Data_n_1_to_300M.zip> <Data_n=300M_plus directory>
"""
import glob
import io
import json
import os
import re
import sys
import zipfile

import numpy as np

Y, GAP = 2_000_000_000, 3000
LINE = re.compile(rb"^\s*(?:n\s*=\s*)?(\d+)", re.M)  # "n | digits | hash" tables and "n = ..." reports


def ints(data):
    return np.array([int(m) for m in LINE.findall(data)], dtype=np.int64)


def main():
    zpath, ddir = sys.argv[1], sys.argv[2]
    bits = np.zeros(Y + 2, dtype=bool)
    files, total_lines = {}, 0
    sources = []
    with zipfile.ZipFile(zpath) as z:
        for e in z.infolist():
            if e.filename.endswith(".txt"):
                sources.append((e.filename, z.read(e)))
    for f in sorted(glob.glob(os.path.join(ddir, "Prime_Q47_*M_*M.txt"))):
        sources.append((os.path.basename(f), open(f, "rb").read()))
    for name, data in sources:
        a = ints(data)
        a = a[(a >= 1) & (a <= Y)]
        dup_before = int(bits[a].sum())
        bits[a] = True
        total_lines += len(a)
        files[name] = {"lines": len(a), "min": int(a.min()) if len(a) else None,
                       "max": int(a.max()) if len(a) else None, "already_seen": dup_before}
        print(name, files[name], flush=True)
        del data, a
    n = np.flatnonzero(bits)
    gaps = np.diff(n)
    big = [(int(n[i]), int(n[i + 1])) for i in np.flatnonzero(gaps > GAP)]
    b = bits[1:Y + 1]
    out = {"files": files, "lines_total": total_lines, "distinct_primes": int(len(n)),
           "doublets": int(np.sum(b[:-1] & b[1:])),
           "triplets": int(np.sum(b[:-2] & b[1:-1] & b[2:])),
           "quadruplets": int(np.sum(b[:-3] & b[1:-2] & b[2:-1] & b[3:])),
           "first_n": int(n[0]), "last_n": int(n[-1]), "gaps_over_%d" % GAP: big}
    print(json.dumps({k: v for k, v in out.items() if k != "files"}, indent=1))
    json.dump(out, open("count_primes_2e9.json", "w"), indent=1)


if __name__ == "__main__":
    main()
