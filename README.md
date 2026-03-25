# Gaia DR3 6656986282721029120 — Dormant Black Hole Candidate

[![arXiv](https://img.shields.io/badge/arXiv-TBD-b31b1b.svg)](https://arxiv.org/abs/TBD)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

Publication-ready analysis of **Gaia DR3 6656986282721029120**, a dormant
black-hole candidate identified via the Gaia DR3 non-single-star (NSS) catalogue.

| Property | Value |
|---|---|
| Gaia DR3 Source ID | 6656986282721029120 |
| RA, Dec (J2016.0) | 288.1005°, −51.8057° |
| *G* magnitude | 6.47 |
| Parallax | 1.470 ± 0.088 mas (d ≈ 680 pc) |
| Orbital period | 443.9 ± 2.5 d |
| Eccentricity | 0.583 ± 0.043 |
| Primary mass (M₁) | ~9.2 M☉ (K-type giant) |
| Companion mass (M₂) | 5.51 M☉ (true mass; Orbital solution) |
| NSS solution type | Orbital (astrometric) |
| Astrometric significance | 28.0σ |

The companion mass of ~5.5 M☉ is well above the neutron-star maximum
(~2.3 M☉) but sits at the lower boundary of the established stellar-mass
black-hole population, making this a **borderline mass-gap / black-hole
candidate**. The primary is a luminous K-type giant, and the companion is
dark — no luminous main-sequence star of compatible mass is detected in the
spectral energy distribution.

## Repository Structure

```
dr3-6656986-black-hole-candidate/
├── scripts/
│   ├── 01_build_target_dataset.py   # Gaia TAP + VizieR + SIMBAD queries
│   ├── 02_fit_sed_extinction.py     # SED fitting, Teff from colour, extinction
│   ├── 03_compute_mass_posterior.py  # Mass posterior & mass-gap analysis
│   ├── 04_companion_exclusion.py    # Luminous-companion exclusion (SED test)
│   ├── 05_alternative_scenarios.py  # WD / NS / MS / stripped-He scenarios
│   ├── 06_make_figures.py           # Publication figures
│   ├── 07_sensitivity_analysis.py   # M₁, parallax, threshold sweeps
│   └── 08_archival_checks.py        # Variability, neighbours, X-ray, kinematics
├── paper/
│   ├── manuscript.tex               # MNRAS-format LaTeX manuscript
│   ├── references.bib               # BibTeX bibliography
│   ├── figures/                     # Generated figures (PDF)
│   └── tables/                      # Generated tables
├── results/                         # JSON output from each script
├── data/
│   └── metadata_notes.md            # Source parameter reference
├── requirements.txt
├── LICENSE
├── CITATION.cff
└── README.md
```

## Quick Start

### Prerequisites

- Python ≥ 3.10
- LaTeX distribution with BibTeX (e.g. MiKTeX, TeX Live)

### Installation

```bash
git clone https://github.com/jayalabaez/dr3-6656986-black-hole-candidate.git
cd dr3-6656986-black-hole-candidate
pip install -r requirements.txt
```

### Run the Analysis

Execute the scripts in order:

```bash
cd scripts
python 01_build_target_dataset.py
python 02_fit_sed_extinction.py
python 03_compute_mass_posterior.py
python 04_companion_exclusion.py
python 05_alternative_scenarios.py
python 06_make_figures.py
python 07_sensitivity_analysis.py
python 08_archival_checks.py
```

### Compile the Manuscript

```bash
cd ../paper
pdflatex manuscript
bibtex manuscript
pdflatex manuscript
pdflatex manuscript
```

## Key Results

- **Companion mass**: M₂ = 5.51 M☉ (true mass from astrometric orbital solution).
- **White dwarf excluded**: 3.8× above Chandrasekhar limit.
- **Neutron star excluded**: 2.4× above TOV maximum.
- **Main-sequence companion disfavoured**: SED composite test finds no blue-band excess
  expected from a ~5.5 M☉ B-type main-sequence star.
- **Mass-gap sensitivity**: M₂ depends on assumed M₁; for M₁ ≥ 7 M☉ the companion
  remains above 3 M☉ (mass-gap object) for all Monte Carlo realisations.
- **High eccentricity** (e = 0.58) is consistent with BH natal kick and argues against
  mass-transfer history.

## Citation

If you use this analysis, please cite:

```bibtex
@ARTICLE{AyalaJoel2025_DR3_6656986,
  author  = {{Ayala}, Joel},
  title   = {{A Borderline Black Hole Candidate from Astrometric Orbital
              Solution with a K-type Giant Primary
              (Gaia DR3 6656986282721029120)}},
  year    = {2025},
  note    = {GitHub repository},
  url     = {https://github.com/jayalabaez/dr3-6656986-black-hole-candidate}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE).

## Acknowledgements

This work has made use of data from the European Space Agency (ESA) mission
*Gaia* (https://www.cosmos.esa.int/gaia), processed by the Gaia Data
Processing and Analysis Consortium (DPAC).
