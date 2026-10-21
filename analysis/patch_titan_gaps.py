#!/usr/bin/env python3
"""
patch_titan_gaps.py -- complete the user's Titan sweeper census of Q53, Q61, Q71 quadruplets (n <= 2e9).

The sweeper scans blocks [s, s + 10^7), s = 1 + k*10^7, but only tests starts i < 10^7 - 3 in each
block, so the three starts s + 10^7 - 3 .. s + 10^7 - 1 per block are never examined.  This script
tests those 600 starts per exponent directly (BPSW on all four values) and re-confirms every reported
quadruplet with BPSW.  Output: titan_census_q53_61_71.json.
"""
import json

import gmpy2

from bh_constants import expected

REPORTED = {
    53: [706668099, 1393928909],
    61: [625149337, 1144583846, 1399766329, 1690426017],
    71: [240500021, 392842707, 682677091, 841911731, 1153177321, 1974811364],
}
Y, CHUNK = 2_000_000_000, 10_000_000


def Q(q, n):
    n = gmpy2.mpz(n)
    return n ** q - (n - 1) ** q


def is_quad(q, n):
    return all(gmpy2.is_bpsw_prp(Q(q, n + j)) for j in range(4))


consts = json.load(open("bh_constants.json"))
out = {}
for q, rep in REPORTED.items():
    gaps = [s + CHUNK - d for s in range(1, Y + 1, CHUNK) for d in (3, 2, 1) if s + CHUNK - d <= Y]
    extra = [n for n in gaps if is_quad(q, n)]
    confirmed = [n for n in rep if is_quad(q, n)]
    quads = sorted(set(confirmed) | set(extra))
    E = expected(consts[str(q)]["4"]["C_partial"]["200000000"], q, 4, Y)
    out[q] = {"reported": rep, "bpsw_confirmed": confirmed, "gap_starts_tested": len(gaps),
              "gap_hits": extra, "quadruplets": quads, "count": len(quads), "expected_BH": E}
    print(q, json.dumps(out[q]), flush=True)
json.dump(out, open("titan_census_q53_61_71.json", "w"), indent=1)
