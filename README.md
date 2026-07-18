# Research Data Supporting "Atomistic mechanisms of hard carbon formation from polyvinylidene chloride"

This repository accompanies the manuscript:

> Litong Wu, Zitong Wu, Zakariya El-Machachi, and Volker L. Deringer,
> "Atomistic mechanisms of hard carbon formation from polyvinylidene chloride" (in preparation).

## Repository contents

The repository is organized into two main directories:

- `data/` stores the training and testing datasets and the final machine-learned interatomic potential.
- `scripts/` contains the analysis and plotting workflows used to generate the manuscript figures.

### `data/`

The current data directory contains:

- The training and testing datasets used to fit and validate the MLIP model:
  - `training-set.xyz`
  - `testing-set.xyz`
- The final trained potential:
  - `pvdc_stagetwo.model`

### `scripts/`

The analysis scripts, processed data, and plotting scripts are organized by manuscript section:

- `1a/` — representative trajectory snapshots at multiple time points, including an OVITO session file for visualization
- `1b/` — temperature-profile analysis
- `1c/` — mass-over-time analysis
- `1d/` — density-over-time analysis
- `2a/` — bond-count analysis
- `2d/` — pre-reactive encounter-time analysis
- `2e/` — event-based lifetime-frequency analysis
- `2f/` — reaction-type analysis
- `3a/` — carbon coordination-cluster analysis
- `3d/` — representative trajectory snapshots at multiple time points with six-membered rings highlighted, together with the corresponding ring-identification script

Each section directory typically contains:

- an analysis script (`*_analysis.py`),
- a plotting script (`*_plot.py`), and
- a generated CSV summary (`*.csv`).

## Data provenance note

A subset of the training and testing datasets were obtained from the following open sources but relabelled:

> **Amorphous carbon and graphitic configurations**: Rowe, P.; Deringer, V. L.; Gasparotto, P.; Csányi, G.; Michaelides, A. *J. Chem. Phys.* **2020**, *153*, 034702.

> **Bulk CH configurations**: Ibragimova, R.; Kuklin, M. S.; Zarrouk, T.; Caro, M. A. *Chem. Mater.* **2025**, *37*, 1094–1110.

> **Hydrocarbon molecules**: Kim, S.; Chen, J.; Cheng, T.; Gindulyte, A.; He, J.; He, S.; Li, Q.; Shoemaker, B. A.;Thiessen, P. A.; Yu, B.; Zaslavsky, L.; Zhang, J.; Bolton, E. E. *Nucleic Acids Res.* **2025**, *53*, D1516–D1525.

> **QM9 molecules**: Wichmann, C.; Ramakrishnan, R.; von Lilienfeld, O. A. *Sci. Data* **2025**, *12*, 9.
