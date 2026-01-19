# Polynomial Prime Quadruplets

**Prime Clustering in Polynomial Q(n) = n^q - (n-1)^q: Quadruplet Distribution and Riemann Zero Correlation**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🔬 Key Discovery

We discovered **15 quadruplet systems** (four consecutive prime-generating integers) for polynomial Q₄₇, where 47 is a **safe prime**. A remarkable correlation (**r = 0.994**) exists between quadruplet positions and Riemann zeta zeros.

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
│   └── plot_figures.py           # Figure generation
├── paper/
│   ├── Q47_paper.tex             # LaTeX source
│   ├── Q47_paper.pdf             # Compiled paper
│   └── Q47_Paper_Final_v5.pdf    # Full paper with figures
└── README.md
```

---

## 🧮 Theoretical Framework

### Subgroup Interference Hypothesis

The interference potential **I(q) = d(q-1) - 2** measures subgroup lattice complexity:

- **Q₄₇:** Linear chain {1} → H₂ → H₂₃ → Z*₄₇ (minimal interference)
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
  title={Prime Clustering in Polynomial Q(n)=n^q-(n-1)^q: Quadruplet Distribution and Riemann Zero Correlation},
  author={Chen, Ruqing},
  year={2026},
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
