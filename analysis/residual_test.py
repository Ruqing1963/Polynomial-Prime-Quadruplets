#!/usr/bin/env python3
"""
residual_test.py -- detrended test of the Q47 quadruplet / Riemann-zero correlation.

Trend removed on both sides:
  quadruplets: a_k = N4(n_k) - (k - 1/2),  N4(Y) = C4 * int_2^Y dt / prod_{i<4} ln Q47(t+i)  (Bateman-Horn)
  zeros:       b_k = Nbar(g_k) - (k - 1/2), Nbar(T) = T/(2 pi) ln(T/(2 pi e)) + 7/8          (Riemann-von Mangoldt)
then Pearson and Spearman correlation of (a_k) with (b_k), with an exact-permutation-style p-value
(2*10^5 random permutations).  Also reproduces the raw r = 0.994 and its sorted-uniform null.
"""
import json
import math

import numpy as np

N = np.array([23159557, 117309848, 136584738, 218787064, 411784485, 423600750, 523331634,
              640399031, 987980498, 1163461515, 1370439187, 1643105964, 1691581855,
              1975860550, 1996430175], float)
G = np.array([14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588,
              37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478,
              52.970321478, 56.446247697, 59.347044003, 60.831778525, 65.112544048])
Q, C4 = 47, json.load(open("bh_constants.json"))["47"]["4"]["C_partial"]["200000000"]


def N4(Y):
    x = np.linspace(math.log(2.0), math.log(Y), 100001)
    f = np.exp(x) / (math.log(Q) + (Q - 1) * x) ** 4
    return C4 * float(np.sum((f[1:] + f[:-1]) / 2 * np.diff(x)))


def Nbar(T):
    return T / (2 * math.pi) * math.log(T / (2 * math.pi * math.e)) + 7 / 8


def spearman(x, y):
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    return np.corrcoef(rx, ry)[0, 1]


k = np.arange(1, 16)
a = np.array([N4(n) for n in N]) - (k - 0.5)
b = np.array([Nbar(g) for g in G]) - (k - 0.5)
r, rho = np.corrcoef(a, b)[0, 1], spearman(a, b)
rng = np.random.default_rng(1)
perm = np.array([np.corrcoef(a, rng.permutation(b))[0, 1] for _ in range(200000)])
p_two = float(np.mean(np.abs(perm) >= abs(r)))

# raw correlation and the sorted-uniform null (best power transform), as in Remark 1
raw = np.corrcoef(N ** (1 / 2.74), G)[0, 1]
alphas = np.linspace(1, 5, 41)
null = []
for _ in range(20000):
    x = np.sort(rng.uniform(1e6, 2e9, 15))
    null.append(max(np.corrcoef(x ** (1 / al), G)[0, 1] for al in alphas))
null = np.array(null)

out = {
    "C4": C4, "N4(2e9)": N4(2e9), "observed": 15,
    "a_k (quadruplet residual)": a.round(3).tolist(),
    "b_k (zero residual)": b.round(3).tolist(),
    "pearson_detrended": r, "spearman_detrended": rho, "perm_p_two_sided": p_two,
    "raw_r_power_2.74": raw, "null_median": float(np.median(null)),
    "null_p95": float(np.percentile(null, 95)), "null_P(r>=raw)": float(np.mean(null >= raw)),
}
print(json.dumps(out, indent=1))
json.dump(out, open("residual_test.json", "w"), indent=1)
