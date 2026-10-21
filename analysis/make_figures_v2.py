#!/usr/bin/env python3
"""Figures for the v2 paper: cumulative Q47 quadruplets vs Bateman-Horn, C4 by exponent, detrended residuals."""
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from bh_constants import expected

BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "legend.frameon": False})
OUT = "../v2/"

# --- Fig: cumulative quadruplet count vs BH ------------------------------------------------
consts = json.load(open("bh_constants.json"))
S = sorted(set(json.load(open("quad_starts_from_zip.json"))) | {
    23159557, 117309848, 136584738, 218787064, 411784485, 423600750, 523331634, 640399031,
    987980498, 1163461515, 1370439187, 1643105964, 1691581855, 1975860550, 1996430175})
S = [n for n in S if n < 8e11]
C4 = consts["47"]["4"]["C_partial"]["200000000"]
grid = np.linspace(2, 8e11, 400)
bh = [expected(C4, 47, 4, y) for y in grid]
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.step(np.array(S) / 1e9, np.arange(1, len(S) + 1), where="post", color=BLUE, lw=2, label="Observed")
ax.plot(grid / 1e9, bh, color=ORANGE, lw=2, ls="--", label="Bateman–Horn, $C_4=6515$")
ax.set_xlabel("$n$ ($\\times 10^9$)")
ax.set_ylabel("Quadruplets with start $\\leq n$")
ax.text(800, len(S), f"  {len(S)}", va="center", color=INK, fontsize=9)
ax.text(800, bh[-1] - 120, f"  {bh[-1]:.0f}", va="center", color=MUTED, fontsize=9)
ax.set_xlim(0, 880)
ax.legend(loc="upper left")
fig.tight_layout()
fig.savefig(OUT + "fig_census_bh.png", dpi=220)

# --- Fig: C4 by exponent, colored by whether 2q+1 is prime -------------------------------
qs = [23, 29, 37, 41, 43, 47, 53, 59, 61, 71]
sg = {23, 29, 41, 53}  # 2q+1 prime
safe = {23, 47, 59}    # (q-1)/2 prime
vals = [consts[str(q)]["4"]["C_partial"]["200000000"] for q in qs]
fig, ax = plt.subplots(figsize=(6.4, 3.4))
x = np.arange(len(qs))
cols = [ORANGE if q in sg else BLUE for q in qs]
ax.bar(x, vals, width=0.62, color=cols, edgecolor="white", linewidth=1.5)
for i, (q, v) in enumerate(zip(qs, vals)):
    ax.text(i, v + 150, f"{v:,.0f}", ha="center", va="bottom", fontsize=8, color=INK)
ax.set_xticks(x, [f"{q}{'*' if q in safe else ''}" for q in qs])
ax.set_xlabel("Exponent $q$   (* = safe prime)")
ax.set_ylabel("Quadruplet constant $C_4(q)$")
ax.grid(axis="x", visible=False)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=BLUE, label="$2q+1$ composite"),
                   Patch(color=ORANGE, label="$2q+1$ prime (Sophie Germain $q$)")], loc="upper left")
ax.set_ylim(0, max(vals) * 1.15)
fig.tight_layout()
fig.savefig(OUT + "fig_c4.png", dpi=220)

# --- Fig: detrended residuals -------------------------------------------------------------
r = json.load(open("residual_test.json"))
a, b = np.array(r["a_k (quadruplet residual)"]), np.array(r["b_k (zero residual)"])
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.9))
k = np.arange(1, 16)
axs[0].plot(k, a, "-o", color=BLUE, lw=2, ms=4, label="Quadruplets $a_k$")
axs[0].plot(k, b, "-s", color=ORANGE, lw=2, ms=4, label="Zeros $b_k$")
axs[0].axhline(0, color=MUTED, lw=0.8)
axs[0].set_xlabel("$k$")
axs[0].set_ylabel("Detrended residual")
axs[0].legend(loc="upper left", fontsize=8)
axs[1].scatter(b, a, s=36, color=BLUE, edgecolor="white", linewidth=1.5, zorder=3)
axs[1].set_xlabel("$b_k$ (zeros)")
axs[1].set_ylabel("$a_k$ (quadruplets)")
axs[1].set_title(f"$r={r['pearson_detrended']:+.3f}$, perm. $p={r['perm_p_two_sided']:.2f}$",
                 fontsize=9, color=INK)
fig.tight_layout()
fig.savefig(OUT + "fig_residuals.png", dpi=220)
print("ok")
