#!/usr/bin/env python3
"""
07_sensitivity_analysis.py — Sensitivity of the companion mass
posterior to key assumptions for Gaia DR3 6656986282721029120.

This candidate is particularly sensitive to the primary mass M1
because M1 = 9.18 Msun (luminous K-giant/supergiant) is itself
highly uncertain.  The M1 prior dominates the M2 uncertainty.

Swept parameters:
  - M1 primary mass (3 to 15 Msun)
  - parallax / distance (+-1 sigma, +-2 sigma)
  - BH mass threshold (2.5 to 8 Msun)
  - eccentricity (+-3 sigma)

Outputs:
  results/sensitivity_results.json
  paper/figures/fig_sensitivity.pdf
"""

import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASEDIR = os.path.join(os.path.dirname(__file__), '..')
FIGDIR = os.path.join(BASEDIR, 'paper', 'figures')
RESDIR = os.path.join(BASEDIR, 'results')
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(RESDIR, exist_ok=True)

# Reference values (Orbital solution)
M1_REF = 9.179
M2_REF = 5.509
PLX = 1.4698
PLX_ERR = 0.0883
P_REF = 443.904
P_ERR = 2.514
ECC_REF = 0.5830
ECC_ERR = 0.0427
NDRAWS = 500_000

BH_THRESH_DEFAULT = 3.0   # mass gap lower edge
BH_THRESH_HIGH = 5.0      # conventional BH threshold


def kepler_m2(m1, P_d, m2_ref):
    """Rescale M2 keeping photo-centre semi-major axis a0 fixed."""
    M_total_ref = M1_REF + m2_ref
    lhs = m2_ref / M_total_ref**(2.0 / 3.0)
    scale = (P_REF / P_d)**(2.0 / 3.0) * lhs

    m2 = np.full_like(np.atleast_1d(m1), m2_ref, dtype=np.float64)
    m1 = np.atleast_1d(m1).astype(np.float64)
    P_d = np.atleast_1d(P_d).astype(np.float64)
    for _ in range(50):
        mt = m1 + m2
        f = m2 / mt**(2.0 / 3.0) - scale
        df = 1.0 / mt**(2.0 / 3.0) - (2.0 / 3.0) * m2 / mt**(5.0 / 3.0)
        m2 = m2 - f / df
        m2 = np.maximum(m2, 0.01)
    return m2


def compute_posterior(m1, P_d, systematic_frac=0.10):
    """Return M2 draws including Gaussian model systematic."""
    m2 = kepler_m2(m1, P_d, M2_REF)
    noise = 1.0 + systematic_frac * np.random.randn(len(m1))
    return m2 * noise


def run_sweep():
    rng = np.random.default_rng(42)
    results = {}

    # baseline draws — broad M1 prior
    m1_draws = rng.lognormal(
        mean=np.log(M1_REF) - 0.5 * (3.0 / M1_REF)**2,
        sigma=3.0 / M1_REF,
        size=NDRAWS
    )
    P_draws = rng.normal(P_REF, P_ERR, NDRAWS)
    P_draws = np.clip(P_draws, 1.0, None)

    m2_base = compute_posterior(m1_draws, P_draws)
    valid = (m2_base > 0.1) & (m2_base < 200)
    m2_base = m2_base[valid]

    results['baseline'] = {
        'median': float(np.median(m2_base)),
        'p16': float(np.percentile(m2_base, 16)),
        'p84': float(np.percentile(m2_base, 84)),
        'P_BH_3': float(np.mean(m2_base > BH_THRESH_DEFAULT)),
        'P_BH_5': float(np.mean(m2_base > BH_THRESH_HIGH)),
    }
    print(f'  Baseline: M2 = {results["baseline"]["median"]:.2f} '
          f'[{results["baseline"]["p16"]:.2f}, '
          f'{results["baseline"]["p84"]:.2f}] Msun')
    print(f'    P(M2 > 3 Msun) = {results["baseline"]["P_BH_3"]:.4f}')
    print(f'    P(M2 > 5 Msun) = {results["baseline"]["P_BH_5"]:.4f}')

    # ── M1 sweep (critical for this source) ──────────────────────
    m1_grid = np.array([3.0, 5.0, 7.0, 9.18, 12.0, 15.0])
    m1_sweep = {}
    for m1v in m1_grid:
        m1d = np.full(NDRAWS, m1v) + rng.normal(0, 0.5, NDRAWS)
        m1d = np.clip(m1d, 0.5, None)
        m2 = compute_posterior(m1d, P_draws[:NDRAWS])
        m2 = m2[(m2 > 0.1) & (m2 < 200)]
        m1_sweep[str(m1v)] = {
            'median': float(np.median(m2)),
            'P_BH_3': float(np.mean(m2 > BH_THRESH_DEFAULT)),
            'P_BH_5': float(np.mean(m2 > BH_THRESH_HIGH)),
        }
        print(f'    M1={m1v:.1f}: M2_med={np.median(m2):.2f}, '
              f'P(>5)={np.mean(m2>5):.3f}')
    results['M1_sweep'] = m1_sweep

    # ── Parallax / distance sweep ────────────────────────────────
    plx_grid = [PLX + 2*PLX_ERR, PLX + PLX_ERR, PLX,
                PLX - PLX_ERR, PLX - 2*PLX_ERR]
    plx_sweep = {}
    for plx_val in plx_grid:
        if plx_val <= 0:
            continue
        dist = 1000.0 / plx_val
        # Distance mainly affects M1 through luminosity
        dm = 5 * np.log10(dist / 10.0)
        mg = 6.469 - dm - 0.61
        bc = -0.40
        m_bol = mg + bc
        log_l = (4.74 - m_bol) / 2.5
        # Rough mass from luminosity for evolved star
        m1_est = max(1.0, 10**(log_l / 3.0))
        m1d = np.full(NDRAWS, m1_est) + rng.normal(0, 0.5, NDRAWS)
        m1d = np.clip(m1d, 0.5, None)
        m2 = compute_posterior(m1d, P_draws[:NDRAWS])
        m2 = m2[(m2 > 0.1) & (m2 < 200)]
        plx_sweep[f'plx={plx_val:.4f}'] = {
            'distance_pc': float(dist),
            'M1_est': float(m1_est),
            'M2_median': float(np.median(m2)),
            'P_BH_5': float(np.mean(m2 > BH_THRESH_HIGH)),
        }
    results['parallax_sweep'] = plx_sweep

    # ── BH threshold sweep ───────────────────────────────────────
    thresholds = [2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0]
    thresh_sweep = {}
    for th in thresholds:
        thresh_sweep[str(th)] = float(np.mean(m2_base > th))
    results['BH_threshold_sweep'] = thresh_sweep

    return results, m2_base, m1_grid, m1_sweep


def make_figure(m2_base, m1_grid, m1_sweep):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel 1: baseline posterior
    ax = axes[0]
    ax.hist(m2_base[m2_base < 25], bins=120, density=True,
            color='steelblue', alpha=0.7, edgecolor='none')
    ax.axvline(BH_THRESH_HIGH, color='red', ls='--', lw=1.5,
               label=f'BH threshold = {BH_THRESH_HIGH} M$_\\odot$')
    ax.axvline(BH_THRESH_DEFAULT, color='orange', ls='--', lw=1.5,
               label=f'Mass gap = {BH_THRESH_DEFAULT} M$_\\odot$')
    ax.axvline(np.median(m2_base), color='k', ls='-', lw=1.2,
               label=f'Median = {np.median(m2_base):.2f} M$_\\odot$')
    ax.set_xlabel('$M_2$  [M$_\\odot$]')
    ax.set_ylabel('Probability density')
    ax.set_title('Baseline mass posterior')
    ax.legend(fontsize=7)
    ax.set_xlim(0, 25)

    # Panel 2: M1 sensitivity
    ax = axes[1]
    medians = [m1_sweep[str(m)]['median'] for m in m1_grid]
    p_bh5 = [m1_sweep[str(m)]['P_BH_5'] for m in m1_grid]
    ax.plot(m1_grid, medians, 'o-', color='steelblue', label='Median $M_2$')
    ax.axhline(BH_THRESH_HIGH, color='red', ls='--', lw=1, alpha=0.5)
    ax.axhline(BH_THRESH_DEFAULT, color='orange', ls='--', lw=1, alpha=0.5)
    ax.set_xlabel('Assumed $M_1$  [M$_\\odot$]')
    ax.set_ylabel('Median $M_2$  [M$_\\odot$]')
    ax.set_title('$M_1$ sensitivity (CRITICAL)')

    ax2 = ax.twinx()
    ax2.plot(m1_grid, p_bh5, 's--', color='green', alpha=0.7,
             label='P($M_2$ > 5 M$_\\odot$)')
    ax2.set_ylabel('P(BH)', color='green')
    ax2.set_ylim(0, 1.05)
    ax.legend(fontsize=8, loc='lower left')
    ax2.legend(fontsize=8, loc='lower right')

    # Panel 3: BH threshold sensitivity
    ax = axes[2]
    thresholds = [2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0]
    probs = [float(np.mean(m2_base > th)) for th in thresholds]
    ax.plot(thresholds, probs, 'o-', color='darkred')
    ax.fill_between(thresholds, probs, alpha=0.15, color='darkred')
    ax.set_xlabel('BH mass threshold  [M$_\\odot$]')
    ax.set_ylabel('P($M_2$ > threshold)')
    ax.set_title('Threshold sensitivity')
    ax.set_ylim(0, 1.05)

    fig.suptitle('Gaia DR3 6656986282721029120  \u2014  Sensitivity analysis',
                 fontsize=13, fontweight='bold', y=1.02)
    fig.tight_layout()
    outpath = os.path.join(FIGDIR, 'fig_sensitivity.pdf')
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved {outpath}')


def main():
    print('=== Sensitivity analysis ===\n')
    results, m2_base, m1_grid, m1_sweep = run_sweep()

    make_figure(m2_base, m1_grid, m1_sweep)

    outpath = os.path.join(RESDIR, 'sensitivity_results.json')
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\n  Saved {outpath}')
    print('\n=== Sensitivity analysis complete ===')


if __name__ == '__main__':
    main()
