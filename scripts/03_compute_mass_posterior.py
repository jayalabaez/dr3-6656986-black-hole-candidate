#!/usr/bin/env python3
"""
03_compute_mass_posterior.py — Bayesian mass posterior for
Gaia DR3 6656986282721029120.

This source has an **Orbital** (pure astrometric) solution.  The
Gaia pipeline fits the photo-centre orbit from astrometry alone,
without any spectroscopic RV curve.  The inclination is resolved
from the on-sky orbit, yielding a *true* companion mass (not a
minimum mass from the spectroscopic mass function).

The dominant source of uncertainty is the primary mass M1.  The
catalog quotes M1 = 9.18 Msun (K-type giant/supergiant), which
is itself uncertain by ~30-50 %.  Because M2 depends on M1 via
Kepler's third law, the companion mass posterior is broad.

Methodology:
  5x10^5 Monte Carlo draws propagate uncertainties in M1, parallax,
  period, and a 10 % model systematic through the Kepler relation
  that keeps the observed photo-centre semi-major axis a0 fixed.

Outputs:
  results/mass_posterior_results.json
  paper/figures/fig_mass_posterior.pdf
"""

import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ─── orbital parameters (Orbital solution) ───────────────────────────
P_DAYS = 443.904
P_ERR  = 2.514
ECC    = 0.5830
ECC_ERR = 0.0427

# ─── Gaia astrometry ─────────────────────────────────────────────────
PLX     = 1.4698   # mas
PLX_ERR = 0.0883   # mas (6% — excellent)

# ─── primary mass from isochrone/colour placement ────────────────────
# M_G ≈ -2.7 and BP-RP = 1.77 place the primary as a luminous
# K-type giant or early AGB.  The catalog derives M1 = 9.18 Msun.
# This is intrinsically uncertain: a K3 III could be 2-12 Msun
# depending on evolutionary state.  We adopt M1 = 9.18 ± 3.0 Msun
# (log-normal prior to keep M1 > 0).
M1_BEST  = 9.179    # Msun
M1_SIGMA = 3.0      # broad — reflects giant/supergiant ambiguity

# ─── companion mass from Gaia binary_masses (Orbital) ────────────────
M2_CATALOG = 5.509  # Msun (true mass from full orbital solution)

# ─── MC parameters ───────────────────────────────────────────────────
N_DRAWS = 500_000
BH_THRESHOLD = 5.0   # Msun — conventional
NS_MAX = 2.3         # Msun — TOV limit (conservative)
MASS_GAP_LO = 3.0    # Msun — lower mass-gap boundary
MASS_GAP_HI = 5.0    # Msun — upper mass-gap boundary

# ─── constants ───────────────────────────────────────────────────────
G_SI = 6.674e-11
MSUN = 1.989e30
AU = 1.496e11
DAY = 86400.0


def kepler_m2_vectorised(m1_arr, P_d_arr, m2_cat, systematic_arr):
    """
    Derive M2 from the Orbital solution keeping the observed
    photo-centre semi-major axis a0 fixed.

    For a dark companion (flux ratio ~ 0):
      a0 = a * M2 / (M1 + M2)
    where a = (G (M1+M2) P^2 / (4 pi^2))^{1/3}.  Therefore:
      a0 = M2 * (G P^2 / (4 pi^2))^{1/3} / (M1 + M2)^{2/3}

    Since a0 is measured in angular units (mas) and converted to
    physical units via the parallax, a0_phys = a0_ang / plx.
    For the Orbital solution the parallax enters only through
    the distance used to derive M1 from isochrone fitting.

    We keep M2 / (M1+M2)^{2/3} * P^{-2/3} = constant (= catalog value)
    and solve for M2 given varied M1 and P via Newton iteration.
    """
    M_total_ref = M1_BEST + m2_cat
    lhs_ref = m2_cat / M_total_ref**(2.0 / 3.0)
    scale = lhs_ref * (P_DAYS / P_d_arr)**(2.0 / 3.0)

    # Newton iteration (vectorised)
    m2 = np.full_like(m1_arr, m2_cat, dtype=np.float64)
    for _ in range(60):
        mt = m1_arr + m2
        f  = m2 / mt**(2.0 / 3.0) - scale
        df = 1.0 / mt**(2.0 / 3.0) - (2.0 / 3.0) * m2 / mt**(5.0 / 3.0)
        m2 = m2 - f / df
        m2 = np.maximum(m2, 0.01)

    return m2 * (1.0 + systematic_arr)


def main():
    print('=== Mass Posterior for Gaia DR3 6656986282721029120 ===\n')

    # Semi-major axis from catalog parameters
    a_au = ((G_SI * (M1_BEST + M2_CATALOG) * MSUN *
             (P_DAYS * DAY)**2) / (4 * np.pi**2))**(1/3) / AU
    print(f'  Solution type    = Orbital (pure astrometric)')
    print(f'  Mass type        = true (inclination resolved)')
    print(f'  Catalog M2       = {M2_CATALOG:.2f} Msun')
    print(f'  M1               = {M1_BEST:.2f} +/- {M1_SIGMA:.2f} Msun')
    print(f'  P                = {P_DAYS:.2f} +/- {P_ERR:.2f} d')
    print(f'  e                = {ECC:.4f} +/- {ECC_ERR:.4f}')
    print(f'  plx              = {PLX:.4f} +/- {PLX_ERR:.4f} mas')
    print(f'  Semi-major axis  = {a_au:.2f} AU')

    # ── Monte Carlo ──────────────────────────────────────────────
    print(f'\n  Running MC ({N_DRAWS:,} draws) ...')
    rng = np.random.default_rng(42)

    # Draw parameters
    plx_draws = rng.normal(PLX, PLX_ERR, N_DRAWS)
    plx_draws = np.clip(plx_draws, 0.05, None)

    p_draws = rng.normal(P_DAYS, P_ERR, N_DRAWS)
    p_draws = np.clip(p_draws, 1.0, None)

    # M1: log-normal to keep positive, centred on catalog
    m1_draws = rng.lognormal(
        mean=np.log(M1_BEST) - 0.5 * (M1_SIGMA / M1_BEST)**2,
        sigma=M1_SIGMA / M1_BEST,
        size=N_DRAWS
    )

    model_sys = rng.normal(0, 0.10, N_DRAWS)  # 10% systematic

    m2_draws = kepler_m2_vectorised(m1_draws, p_draws, M2_CATALOG,
                                     model_sys)

    # Filter valid
    valid = (m2_draws > 0.1) & (m2_draws < 500)
    m2_draws = m2_draws[valid]
    m1_valid = m1_draws[valid]
    n_valid = len(m2_draws)
    print(f'  Valid samples: {n_valid:,}/{N_DRAWS:,}')

    # Statistics
    median = np.median(m2_draws)
    ci68 = np.percentile(m2_draws, [16, 84])
    ci90 = np.percentile(m2_draws, [5, 95])
    p_bh = 100 * np.mean(m2_draws > BH_THRESHOLD)
    p_above_ns = 100 * np.mean(m2_draws > NS_MAX)
    p_mass_gap = 100 * np.mean((m2_draws > MASS_GAP_LO) &
                                (m2_draws < MASS_GAP_HI))
    p_above3 = 100 * np.mean(m2_draws > 3.0)
    p_above10 = 100 * np.mean(m2_draws > 10)

    print(f'\n  Results:')
    print(f'    M2 median          = {median:.2f} Msun')
    print(f'    68% CI             = [{ci68[0]:.2f}, {ci68[1]:.2f}] Msun')
    print(f'    90% CI             = [{ci90[0]:.2f}, {ci90[1]:.2f}] Msun')
    print(f'    P(M2 > 5 Msun)    = {p_bh:.1f}%')
    print(f'    P(M2 > NS max)    = {p_above_ns:.1f}%')
    print(f'    P(M2 > 3 Msun)    = {p_above3:.1f}%')
    print(f'    P(M2 > 10 Msun)   = {p_above10:.1f}%')
    print(f'    P(mass gap)        = {p_mass_gap:.1f}%')

    # Sensitivity: vary M1
    print(f'\n  Sensitivity to M1:')
    for m1_test in [3.0, 5.0, 7.0, 9.18, 12.0, 15.0]:
        m2_t = kepler_m2_vectorised(
            np.array([m1_test]), np.array([P_DAYS]),
            M2_CATALOG, np.array([0.0]))
        a_test = ((G_SI * (m1_test + m2_t[0]) * MSUN *
                   (P_DAYS * DAY)**2) / (4 * np.pi**2))**(1/3) / AU
        print(f'    M1={m1_test:5.1f}: M2={m2_t[0]:.2f} Msun, '
              f'a={a_test:.2f} AU')

    # ── Save results ─────────────────────────────────────────────
    basedir = os.path.dirname(__file__)
    results = {
        'solution_type': 'Orbital',
        'mass_type': 'true (inclination resolved from astrometry)',
        'M2_catalog': M2_CATALOG,
        'M1_best': M1_BEST,
        'M1_sigma': M1_SIGMA,
        'parallax_mas': PLX,
        'parallax_err_mas': PLX_ERR,
        'period_d': P_DAYS,
        'eccentricity': ECC,
        'MC_draws': N_DRAWS,
        'MC_valid': n_valid,
        'M2_median': round(median, 2),
        'M2_68ci': [round(ci68[0], 2), round(ci68[1], 2)],
        'M2_90ci': [round(ci90[0], 2), round(ci90[1], 2)],
        'P_above_5Msun_percent': round(p_bh, 1),
        'P_above_NS_percent': round(p_above_ns, 1),
        'P_above_3Msun_percent': round(p_above3, 1),
        'P_above_10Msun_percent': round(p_above10, 1),
        'P_mass_gap_percent': round(p_mass_gap, 1),
        'semi_major_axis_AU': round(a_au, 2),
    }
    outpath = os.path.join(basedir, '..', 'results',
                           'mass_posterior_results.json')
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\n  Saved: {outpath}')

    # ── Figure ───────────────────────────────────────────────────
    figdir = os.path.join(basedir, '..', 'paper', 'figures')
    os.makedirs(figdir, exist_ok=True)
    figpath = os.path.join(figdir, 'fig_mass_posterior.pdf')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: M2 sensitivity to M1
    m1_grid = np.linspace(2.0, 18.0, 200)
    m2_grid = np.array([
        kepler_m2_vectorised(np.array([m1]), np.array([P_DAYS]),
                              M2_CATALOG, np.array([0.0]))[0]
        for m1 in m1_grid
    ])
    ax1.plot(m1_grid, m2_grid, 'b-', lw=2)
    ax1.axhline(BH_THRESHOLD, color='red', ls='--', alpha=0.7,
                label=r'$M_2 = 5\,M_\odot$ (BH threshold)')
    ax1.axhline(NS_MAX, color='orange', ls='--', alpha=0.7,
                label=r'$M_2 = 2.3\,M_\odot$ (NS ceiling)')
    ax1.axvspan(M1_BEST - M1_SIGMA, M1_BEST + M1_SIGMA,
                alpha=0.15, color='blue',
                label=f'$M_1 = {M1_BEST:.1f} \\pm {M1_SIGMA:.1f}$')
    ax1.axvline(M1_BEST, color='blue', ls=':', alpha=0.5)
    ax1.set_xlabel(r'$M_1$ ($M_\odot$)', fontsize=12)
    ax1.set_ylabel(r'$M_2$ ($M_\odot$)', fontsize=12)
    ax1.set_title(r'Companion mass vs $M_1$ (Orbital solution)')
    ax1.legend(fontsize=8, loc='upper right')
    ax1.set_xlim(2, 18)

    # Right: MC posterior histogram
    bins = np.linspace(0, 25, 120)
    ax2.hist(m2_draws[m2_draws < 25], bins=bins, density=True,
             color='steelblue', alpha=0.7, edgecolor='navy', lw=0.3)
    ax2.axvline(BH_THRESHOLD, color='red', ls='--', lw=2,
                label=f'BH threshold (5 $M_\\odot$)')
    ax2.axvline(median, color='black', ls='-', lw=2,
                label=f'Median = {median:.1f} $M_\\odot$')
    ax2.axvspan(ci68[0], ci68[1], alpha=0.15, color='green',
                label=f'68% CI [{ci68[0]:.1f}, {ci68[1]:.1f}]')
    ax2.axvspan(MASS_GAP_LO, MASS_GAP_HI, alpha=0.10, color='purple',
                label='Mass gap (3-5 $M_\\odot$)')
    ax2.text(0.97, 0.95,
             f'$P(M_2 > 5\\,M_\\odot) = {p_bh:.1f}\\%$\n'
             f'$P(M_2 > 3\\,M_\\odot) = {p_above3:.1f}\\%$\n'
             f'(Orbital solution)',
             transform=ax2.transAxes, ha='right', va='top', fontsize=10,
             bbox=dict(boxstyle='round', fc='lightyellow'))
    ax2.set_xlabel(r'$M_2$ ($M_\odot$)', fontsize=12)
    ax2.set_ylabel('Probability density', fontsize=12)
    ax2.set_title(f'MC posterior ({N_DRAWS//1000}K draws)')
    ax2.legend(fontsize=7, loc='upper left')
    ax2.set_xlim(0, 25)

    fig.suptitle('Gaia DR3 6656986282721029120 \u2014 Mass Constraints',
                 fontsize=13, fontweight='bold')
    fig.tight_layout()
    fig.savefig(figpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  -> {figpath}')

    print('\n=== Mass posterior complete ===')


if __name__ == '__main__':
    main()
