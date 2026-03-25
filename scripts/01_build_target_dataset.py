#!/usr/bin/env python3
"""
01_build_target_dataset.py — Automated data retrieval for
Gaia DR3 6656986282721029120 from Gaia TAP, VizieR, and SIMBAD.

Outputs:
  results/gaia_query_results.json
"""

import json, os, sys

SOURCE_ID = 6656986282721029120
RA = 288.1005
DEC = -51.8057
SEARCH_RADIUS = 5.0  # arcsec

BASEDIR = os.path.join(os.path.dirname(__file__), '..')
RESDIR = os.path.join(BASEDIR, 'results')
os.makedirs(RESDIR, exist_ok=True)

# ── Attempt imports ──────────────────────────────────────────────────
try:
    from astroquery.gaia import Gaia
    from astroquery.vizier import Vizier
    from astroquery.simbad import Simbad
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    HAS_ASTROQUERY = True
except ImportError:
    HAS_ASTROQUERY = False


def query_gaia():
    """Query Gaia DR3 main, NSS orbit, and binary_masses tables."""
    job = Gaia.launch_job(f"""
    SELECT g.source_id, g.ra, g.dec, g.parallax, g.parallax_error,
           g.pmra, g.pmdec, g.phot_g_mean_mag, g.bp_rp,
           g.phot_bp_mean_mag, g.phot_rp_mean_mag,
           g.ruwe, g.astrometric_excess_noise_sig,
           g.radial_velocity, g.radial_velocity_error,
           g.teff_gspphot, g.logg_gspphot
    FROM gaiadr3.gaia_source AS g
    WHERE g.source_id = {SOURCE_ID}
    """)
    main_row = job.get_results()

    job2 = Gaia.launch_job(f"""
    SELECT *
    FROM gaiadr3.nss_two_body_orbit AS n
    WHERE n.source_id = {SOURCE_ID}
    """)
    nss_row = job2.get_results()

    job3 = Gaia.launch_job(f"""
    SELECT *
    FROM gaiadr3.binary_masses AS b
    WHERE b.source_id = {SOURCE_ID}
    """)
    mass_row = job3.get_results()

    return main_row, nss_row, mass_row


def query_vizier():
    """Query 2MASS, AllWISE, GALEX, ROSAT, XMM via VizieR."""
    coord = SkyCoord(ra=RA, dec=DEC, unit='deg', frame='icrs')
    radius = SEARCH_RADIUS * u.arcsec
    results = {}

    for cat, name in [('II/246', '2MASS'), ('II/328', 'AllWISE'),
                       ('II/335', 'GALEX'), ('IX/47', 'ROSAT_2RXS'),
                       ('IX/68', 'XMM_4XMM')]:
        try:
            v = Vizier(columns=['**'], row_limit=5)
            tbl = v.query_region(coord, radius=radius, catalog=cat)
            if tbl:
                results[name] = {c: str(tbl[0][c][0]) for c in tbl[0].colnames[:20]}
            else:
                results[name] = 'No match'
        except Exception as e:
            results[name] = f'Query failed: {e}'

    return results


def query_simbad():
    """Query SIMBAD for object type and cross-IDs."""
    coord = SkyCoord(ra=RA, dec=DEC, unit='deg', frame='icrs')
    custom = Simbad()
    custom.add_votable_fields('otype', 'sp', 'flux(V)', 'ids')
    try:
        result = custom.query_region(coord, radius=5 * u.arcsec)
        if result:
            return {c: str(result[c][0]) for c in result.colnames[:15]}
    except Exception as e:
        return f'SIMBAD query failed: {e}'
    return 'No SIMBAD match'


def fallback_summary():
    """Return catalog parameters when online queries are unavailable."""
    return {
        'source_id': SOURCE_ID,
        'ra': RA,
        'dec': DEC,
        'l': 345.358,
        'b': -24.119,
        'parallax': 1.4698,
        'parallax_error': 0.0883,
        'distance_pc': 680.4,
        'G': 6.469,
        'BP': 7.329,
        'RP': 5.557,
        'BP_RP': 1.773,
        'teff_gspphot': None,
        'logg_gspphot': None,
        'rv_km_s': -9.825,
        'rv_err': 2.500,
        'ruwe': 5.376,
        'en_sig': 691.7,
        'period_d': 443.904,
        'period_err': 2.514,
        'eccentricity': 0.5830,
        'ecc_err': 0.0427,
        'significance': 27.98,
        'gof': 8.17,
        'M1_msun': 9.179,
        'M2_msun': 5.509,
        'sol_type': 'Orbital',
        'mass_type': 'true',
        'score': 80,
        'category': 'BH',
    }


def main():
    print('=== Building target dataset ===\n')
    print(f'  Source: Gaia DR3 {SOURCE_ID}')
    print(f'  RA, Dec = {RA:.4f}, {DEC:.4f}\n')

    data = {'target': f'Gaia DR3 {SOURCE_ID}'}

    if HAS_ASTROQUERY:
        try:
            print('  Querying Gaia TAP ...')
            main_row, nss_row, mass_row = query_gaia()
            data['gaia_main'] = {c: str(main_row[c][0]) for c in main_row.colnames}
            if len(nss_row) > 0:
                data['gaia_nss'] = {c: str(nss_row[c][0]) for c in nss_row.colnames}
            if len(mass_row) > 0:
                data['gaia_masses'] = {c: str(mass_row[c][0]) for c in mass_row.colnames}
            print('  Querying VizieR ...')
            data['vizier'] = query_vizier()
            print('  Querying SIMBAD ...')
            data['simbad'] = query_simbad()
        except Exception as e:
            print(f'  Online queries failed: {e}')
            data['fallback'] = fallback_summary()
    else:
        print('  astroquery not available — using catalog fallback')
        data['fallback'] = fallback_summary()

    outpath = os.path.join(RESDIR, 'gaia_query_results.json')
    with open(outpath, 'w') as f:
        json.dump(data, f, indent=2, default=str)
    print(f'\n  Saved: {outpath}')
    print('\n=== Dataset build complete ===')


if __name__ == '__main__':
    main()
