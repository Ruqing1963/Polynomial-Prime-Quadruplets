# Polynomial Prime Quadruplets

**Prime Clustering in Polynomial Q(n) = n^q - (n-1)^q: Quadruplet Distribution and Riemann Zero Correlation**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23517573-blue.svg)](https://doi.org/10.5281/zenodo.23517573)

## 🆕 Version 2 (October 2026)

**A 200,000-digit Q47 probable prime:**

```
Q47(586135293385*10^4336 + 1647) = (586135293385*10^4336+1647)^47 - (586135293385*10^4336+1646)^47
```

* Exactly 200,000 digits (664,383 bits); 3-, 5-, 7- and 11-PRP (PFGW 4.1.7 / GWNUM 31.4)
* Found by a cyclotomic sieve (prime factors ≡ 1 mod 94, 46 classes per prime) over 150,000 values of n, sieved to 2×10^10
* The discovery rate matches Bateman–Horn. This is a PRP, not a proven prime.
* Details, logs and the search script are in [`prp200k/`](prp200k/). The paper is [`paper/Q47_paper_v2.pdf`](paper/Q47_paper_v2.pdf), Section 7.

Other changes in v2:
* **Correction:** the subgroup lattice of Z*₄₇ is a diamond {1} ⊂ H₂, H₂₃ ⊂ Z*₄₇, not a chain (H₂ ⊄ H₂₃ since 2 ∤ 23). I(47) = 2 is unchanged.
* **Calibration:** against a null model of 15 randomly placed, sorted positions, r ≥ 0.994 occurs about 4.7% of the time (median 0.985), so the Riemann-zero correlation is marginal (P ≈ 0.05). See Remark 1 of v2.

## 🔬 Key Discovery

We discovered **15 quadruplet systems** (four consecutive prime-generating integers) for polynomial Q₄₇, where 47 is a **safe prime**. A correlation (**r = 0.994**) exists between quadruplet positions and Riemann zeta zeros (see the v2 null-model calibration above).

### Core Formula

```
n_k^(1/2.74) ∝ γ_k    where r = 0.994
```

### Effective Modulus

```
q_eff = e^2.74 = 15.5 ± 0.5 ≈ q/3 = 15.67
```

This suggests a **triplet-induced dimensional reduction** in the sieve space.

---

## 📊 Q47 Quadruplet Coordinates (Verified)

| k  | Starting Position n | Range  |
|----|---------------------|--------|
| 1  | 23,159,557          | 0.02B  |
| 2  | 117,309,848         | 0.12B  |
| 3  | 136,584,738         | 0.14B  |
| 4  | 218,787,064         | 0.22B  |
| 5  | 411,784,485         | 0.41B  |
| 6  | 423,600,750         | 0.42B  |
| 7  | 523,331,634         | 0.52B  |
| 8  | 640,399,031         | 0.64B  |
| 9  | 987,980,498         | 0.99B  |
| 10 | 1,163,461,515       | 1.16B  |
| 11 | 1,370,439,187       | 1.37B  |
| 12 | 1,643,105,964       | 1.64B  |
| 13 | 1,691,581,855       | 1.69B  |
| 14 | 1,975,860,550       | 1.98B  |
| 15 | 1,996,430,175       | 2.00B  |

---

## 📈 Seven-Polynomial Comparison

| Q_q | Safe Prime | I(q) | S₃      | Quadruplets (n<2B) |
|-----|------------|------|---------|-------------------|
| Q₃₇ | No         | 7    | +20.4%  | 2 (decay)         |
| Q₄₁ | No         | 6    | +48.8%  | 0                 |
| Q₄₃ | No         | 6    | +15.3%  | 1 (decay)         |
| **Q₄₇** | **Yes** | **2** | +35.0% | **15 (growth)** |
| Q₅₃ | No         | 4    | +30.7%  | 0                 |
| Q₆₁ | No         | 10   | -16.7%  | 0                 |
| Q₇₁ | No         | 6    | +20.5%  | 0                 |

**Key Finding:** Q₄₇ is the ONLY polynomial showing quadruplet growth with extended range.

---

## 📁 Repository Structure

```
├── data/
│   ├── quadruplets_Q47.json      # Complete quadruplet coordinates
│   ├── quadruplets_Q47.csv       # CSV format with Riemann comparison
│   ├── seven_polynomials.json    # Seven-polynomial statistics
│   └── riemann_zeros_first15.csv # Riemann zeta zeros
├── figures/
│   ├── fig1_riemann_correlation.pdf
│   ├── fig2_seven_polynomials.pdf
│   ├── fig3_quadruplet_positions.pdf
│   └── fig4_lattice_comparison.pdf
├── scripts/
│   ├── analyze_quadruplets.py    # Correlation analysis
│   ├── plot_figures.py           # Figure generation
│   └── plot_lattice_v2.py        # Corrected subgroup-lattice figure (v2)
├── paper/
│   ├── Q47_paper.tex             # v1 LaTeX source
│   ├── Q47_paper.pdf             # v1 compiled paper
│   ├── Q47_paper_v2.tex          # v2 LaTeX source
│   ├── Q47_paper_v2.pdf          # v2 compiled paper
│   └── fig*.png                  # v2 figures
├── prp200k/                      # 200,000-digit PRP: search script, sieve state, PFGW logs,
│                                 #   remaining/ = untested candidates + PFGW checkpoints
└── README.md
```

---

## 🧮 Theoretical Framework

### Subgroup Interference Hypothesis

The interference potential **I(q) = d(q-1) - 2** measures subgroup lattice complexity:

- **Q₄₇:** Minimal diamond {1} ⊂ H₂, H₂₃ ⊂ Z*₄₇ (minimal interference)
- **Q₄₁:** Branching network (high interference)

### Asymptotic Conjecture

The relative quadruplet density **R(n) → ∞** as n → ∞:

```
R(n) = P_Q47(n) / P_rnd(n) ~ (ln n)^Δ → ∞
```

---

## 📚 Citation

```bibtex
@article{chen2026polynomial,
  title={Prime Clustering in Polynomial Q(n)=n^q-(n-1)^q: Quadruplet Distribution, Riemann Zero Correlation, and a 200,000-Digit Q47 Probable Prime},
  author={Chen, Ruqing},
  year={2026},
  note={Version 2},
  doi={10.5281/zenodo.23517573},
  institution={GUT Geoservice Inc., Montreal}
}
```

---

## 📧 Contact

**Ruqing Chen**  
GUT Geoservice Inc., Montreal, Canada  
Email: ruqing@hotmail.com

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🔗 References

1. Bateman & Horn (1962). *Math. Comp.* 16(79), 363-367.
2. Hardy & Littlewood (1923). *Acta Math.* 44, 1-70.
3. Montgomery (1973). *Proc. Symp. Pure Math.* 24, 181-193.
4. Odlyzko (1987). *Math. Comp.* 48(177), 273-308.
