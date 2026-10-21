# Polynomial Prime Quadruplets

**Prime Constellations of Q(n) = n^q − (n−1)^q: A Bateman–Horn Census to 8×10¹¹, Certified 1000-Digit Quadruplets, and a 200,000-Digit Q47 Probable Prime**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23517573-blue.svg)](https://doi.org/10.5281/zenodo.23517573)

Paper: [`paper/Q47_paper_v2.pdf`](paper/Q47_paper_v2.pdf) (Version 2, October 2026).

## Summary of Version 2

Q_q(n) = n^q − (n−1)^q is the degree-(q−1) analogue of the cuban primes (q = 3). Every prime factor of Q_q(n) is ≡ 1 (mod 2q). Version 2 tests the Bateman–Horn conjecture on the prime values and prime k-tuples Q(n), Q(n+1), …, Q(n+k−1). It finds agreement at every scale tested.

| Q47 data set | Observed | Bateman–Horn |
|---|---|---|
| Primes, n ≤ 2×10⁹ | 18,473,624 | 18,464,550 |
| Doublets / triplets, n ≤ 2×10⁹ | 176,896 / 1,779 | 176,780 / 1,752 |
| Quadruplets, n < 8×10¹¹ (each BPSW-verified) | 2,393 | 2,404.3 |
| Quintuplets, n < 8×10¹¹ | 21 | 17.9 |
| Sextuplets, n < 8×10¹¹ | 0 | 0.14 |
| Primes in 10⁴⁴ < n ≤ 10⁴⁴+10⁶ (2026 digits; independently rescanned) | 1,870 | 1,861.6 |
| Primes in 10⁷ consecutive n near 1.004×10⁴⁴ (independently rescanned) | 18,715 | 18,614 |

* **The C₄ constants explain the differences between exponents.** The constants C_k(q) were computed for ten exponents, k ≤ 6. The quadruplet constant C₄(q) depends mainly on the smallest prime p ≡ 1 (mod q). When 2q+1 is prime (q = 23, 29, 41, 53), C₄ drops by an order of magnitude. Whether q is a safe prime plays no role.
* **Correction of Version 1.** A complete recount of quadruplets with n ≤ 2×10⁹ for seven exponents finds
  * 37: 23
  * 41: 3
  * 43: 10
  * 47: 15
  * 53: 2
  * 61: 4
  * 71: 6

  Version 1 reported 2, 0, 1, 15, 0, 0, 0. Every recount agrees with Bateman–Horn (|z| < 1.2). The apparent Q47 anomaly in Version 1 came from incomplete counts for the other exponents.
* **Control search.** For eight exponents (n ≤ 10⁷) the counts of primes, doublets and triplets match Bateman–Horn. The two safe primes (23, 59) show no surplus.
* **Riemann zeros.** After detrending both sequences, the correlation between quadruplet positions and Riemann zeros is r = −0.0003 (permutation p = 0.999). The earlier r = 0.994 came from the shared growth trend. The q_eff = q/3 rule is withdrawn.
* **1000-digit quadruplets.** There are 11 quadruplets near n ≈ 4.8×10²¹. Eight of them have all four primes certified by ECPP (Primo); see [`certificates/`](certificates/).
* **200,000-digit probable prime.**

  ```
  P* = Q47(586135293385*10^4336 + 1647)
  ```

  It is exactly 200,000 digits long. It passes Fermat tests to bases 3, 5, 7 and 11 (PFGW), and strong base-11 and BPSW tests (GMP). It is a probable prime, not a proven prime. See [`prp200k/`](prp200k/).
* **Lattice correction.** The subgroup lattice of (Z/47Z)* is a diamond, not a chain.

## Q47 quadruplets with n ≤ 2×10⁹

| k | n | k | n | k | n |
|---|---|---|---|---|---|
| 1 | 23,159,557 | 6 | 423,600,750 | 11 | 1,370,439,187 |
| 2 | 117,309,848 | 7 | 523,331,634 | 12 | 1,643,105,964 |
| 3 | 136,584,738 | 8 | 640,399,031 | 13 | 1,691,581,855 |
| 4 | 218,787,064 | 9 | 987,980,498 | 14 | 1,975,860,550 |
| 5 | 411,784,485 | 10 | 1,163,461,515 | 15 | 1,996,430,175 |

All 2,393 starts below 8×10¹¹ are in [`data/quadruplets_Q47_lt_8e11.csv`](data/quadruplets_Q47_lt_8e11.csv) (also `.json`).

## Repository structure

```
├── paper/
│   ├── Q47_paper_v2.tex / .pdf     # Version 2 (current)
│   ├── Q47_paper.tex / .pdf        # Version 1 (superseded)
│   └── fig*.png                    # figures
├── analysis/                       # Version 2 programs and results (Python 3, numpy, gmpy2)
│   ├── bh_constants.py             # Bateman–Horn constants C_k(q) -> bh_constants*.json
│   ├── quad_census.py              # segmented-sieve quadruplet census (q = 37, 41, 43, 47)
│   ├── patch_titan_gaps.py         # block-boundary patch + BPSW for q = 53, 61, 71
│   ├── control_search.py           # primes / k-tuples for eight exponents, n ≤ 10^7
│   ├── count_primes_2e9.py         # Q47 prime count n ≤ 2e9 from the search files
│   ├── patch_prime_gaps.py         # re-test of data gaps (adds 53 missed primes)
│   ├── rescan_window.py            # independent rescan of the 10^44 windows (sieve + PFGW + BPSW)
│   ├── verify_all_quads.py         # BPSW check of all 2,393 quadruplets
│   ├── verify_tuples.py            # quintuplets and large windows
│   ├── residual_test.py            # detrended Riemann-zero test
│   └── make_figures_v2.py          # figures
├── certificates/primo_1000digit/   # Primo ECPP inputs (.in) and certificates (.out)
├── prp200k/                        # 200,000-digit PRP: search, sieve state, logs, remaining candidates
├── data/quadruplets_Q47_lt_8e11.*  # all 2,393 Q47 quadruplet starts n < 8e11
├── data/, figures/, scripts/       # other files: Version 1 material
└── README.md
```

## Citation

```bibtex
@misc{chen2026q47v2,
  title  = {Prime Constellations of $Q(n)=n^q-(n-1)^q$: A Bateman--Horn Census to $8\times10^{11}$,
            Certified 1000-Digit Quadruplets, and a 200,000-Digit $Q_{47}$ Probable Prime},
  author = {Chen, Ruqing},
  year   = {2026},
  note   = {Version 2},
  doi    = {10.5281/zenodo.23517573}
}
```

Earlier work on Q47: Zenodo DOIs 10.5281/zenodo.18526470, 18701355, 18728540, 20753750, 22177596.

## Contact

**Ruqing Chen**, GUT Geoservice Inc., Montreal, Canada. Email: ruqing@hotmail.com

## License

MIT. See [LICENSE](LICENSE).
