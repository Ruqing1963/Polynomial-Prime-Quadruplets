#!/usr/bin/env python3
"""
Generate publication-quality figures for Q47 paper.
Author: Ruqing Chen
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Data
QUADRUPLETS = np.array([
    23159557, 117309848, 136584738, 218787064, 411784485,
    423600750, 523331634, 640399031, 987980498, 1163461515,
    1370439187, 1643105964, 1691581855, 1975860550, 1996430175
])

RIEMANN_ZEROS = np.array([
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918720, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544
])

def plot_riemann_correlation():
    """Figure 1: Riemann zero correlation."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: Power transform
    alpha = 2.74
    n_transformed = QUADRUPLETS ** (1/alpha)
    
    ax1 = axes[0]
    ax1.scatter(n_transformed, RIEMANN_ZEROS, s=100, c='#2E86AB', edgecolors='black', zorder=5)
    
    # Fit line
    slope, intercept, r, _, _ = stats.linregress(n_transformed, RIEMANN_ZEROS)
    x_fit = np.linspace(n_transformed.min(), n_transformed.max(), 100)
    ax1.plot(x_fit, slope * x_fit + intercept, 'r--', lw=2, label=f'r = {r:.3f}')
    
    ax1.set_xlabel(r'$n^{1/2.74}$', fontsize=12)
    ax1.set_ylabel(r'$\gamma_k$ (Riemann zero)', fontsize=12)
    ax1.set_title('Optimal Power Transform', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # Right: Log-linear
    log_n = np.log(QUADRUPLETS)
    
    ax2 = axes[1]
    ax2.scatter(log_n, RIEMANN_ZEROS, s=100, c='#28A745', edgecolors='black', zorder=5)
    
    slope2, intercept2, r2, _, _ = stats.linregress(log_n, RIEMANN_ZEROS)
    x_fit2 = np.linspace(log_n.min(), log_n.max(), 100)
    ax2.plot(x_fit2, slope2 * x_fit2 + intercept2, 'r--', lw=2, label=f'r = {r2:.3f}')
    
    ax2.set_xlabel(r'$\ln(n)$', fontsize=12)
    ax2.set_ylabel(r'$\gamma_k$ (Riemann zero)', fontsize=12)
    ax2.set_title('Log-Linear Fit', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('../figures/fig1_riemann_correlation.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('../figures/fig1_riemann_correlation.png', dpi=150, bbox_inches='tight')
    print("✓ Figure 1 saved")
    plt.close()


def plot_quadruplet_positions():
    """Figure 3: Quadruplet position distribution."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    k = np.arange(1, 16)
    n_billions = QUADRUPLETS / 1e9
    
    colors = ['#28A745' if i in [1, 2, 3, 13, 14] else '#2E86AB' for i in range(15)]
    
    bars = ax.barh(k, n_billions, color=colors, edgecolor='black', height=0.7)
    
    ax.set_xlabel('Position n (billions)', fontsize=12)
    ax.set_ylabel('Quadruplet index k', fontsize=12)
    ax.set_title('Q47 Quadruplet Distribution (n ≤ 2×10⁹)', fontsize=14, fontweight='bold')
    ax.set_yticks(k)
    ax.invert_yaxis()
    ax.grid(True, axis='x', alpha=0.3)
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#28A745', edgecolor='black', label='Newly discovered'),
        Patch(facecolor='#2E86AB', edgecolor='black', label='Previously known')
    ]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig('../figures/fig3_quadruplet_positions.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('../figures/fig3_quadruplet_positions.png', dpi=150, bbox_inches='tight')
    print("✓ Figure 3 saved")
    plt.close()


if __name__ == "__main__":
    import os
    os.makedirs('../figures', exist_ok=True)
    
    plot_riemann_correlation()
    plot_quadruplet_positions()
    print("\n✓ All figures generated")
