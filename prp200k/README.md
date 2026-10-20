# A 200,000-digit Q47 probable prime

```
P* = Q47(n*) = n*^47 - (n*-1)^47,   n* = 586135293385*10^4336 + 1647
```

| Property | Value |
|---|---|
| Decimal digits | 200,000 exactly |
| Bits | 664,383 |
| Leading digits | `1000000000069212334213...` |
| Trailing digits | `...716306856997037706148011906447` |
| SHA-256 of the decimal string | `bd23db83e0c523fd4ddbf4a0a885bd0986db93eacf6fd2051e4a887d1f0a8c70` |
| P* mod 282 | 1 (as for every value of Q47) |
| Status | **Probable prime (PRP)**, not proven |

## Verification

PFGW 4.1.7 (64-bit, build 20251229), GWNUM 31.4, Windows 11.

| Base | Command | Result | Time |
|---|---|---|---|
| 3 | `pfgw64 -f <file>` | 3-PRP | 277.5 + 89.6 s |
| 5 | `pfgw64 -b5 -f <file>` | 5-PRP | 274.3 + 88.7 s |
| 7 | `pfgw64 -b7 -q"..."` | 7-PRP | 280.3 s |
| 11 | `pfgw64 -b11 -q"..."` | 11-PRP | 320.7 s |

PFGW input line:

```
(586135293385*10^4336+1647)^47-(586135293385*10^4336+1647-1)^47
```

Screenshots of the base 3, 5 and 7 runs are in `screenshots/`. `gmp_verify.py` prints the SHA-256 hash and runs an independent strong-PRP test (base 11) with GMP through gmpy2; it takes much longer than PFGW.

## Search

`q47.py` sieves and drives the search. Every prime factor of Q47(n) is 1 mod 94, and each such prime p removes the 46 classes n ≡ 1/(1-u) (mod p), with u running over the nontrivial 47th roots of unity mod p.

* Window: n = 586135293385*10^4336 + t, 0 ≤ t < 150,000
* Sieved to p ≤ 2×10^10, which leaves 30,809 candidates (state file `q47.npz`)
* Bateman–Horn: P(prime | survives sieve) ≈ e^γ ln B / ln Q ≈ 1/10,900, so the window should hold about 2.8 PRPs
* P* (t = 1647) was the 341st survivor of the first of 8 parallel slices. About 2.7×10^3 candidates had been tested at that point, with 0.25 PRPs expected.

`pfgw_logs/` has the PRP hit log and the composite results (with RES64 residues) for each slice.

## Reproduce

```
python q47.py init   --digits 200000 --length 150000 --state q47.npz
python q47.py sieve  --state q47.npz --pmax 2e10
python q47.py export --state q47.npz --split 8 --outdir pfgw_run
python q47.py selftest
```

See Section 7 of `paper/Q47_paper_v2.pdf` (DOI [10.5281/zenodo.23517573](https://doi.org/10.5281/zenodo.23517573)).
