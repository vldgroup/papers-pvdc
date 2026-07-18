# Research Data Supporting "Atomistic mechanisms of hard carbon formation from polyvinylidene chloride"

This repository accompanies the manuscript:

> Litong Wu, Zitong Wu, Zakariya El-Machachi, and Volker L. Deringer,
> "Atomistic mechanisms of hard carbon formation from polyvinylidene chloride" (in preparation).


The repository is organized into two main directories:

- `data/` stores the final machine-learned interatomic potential and the training and testing datasets.
- `scripts/` contains the analysis and plotting scripts used to generate the manuscript figures.

### `data/`

The data directory contains:

- The final trained potential:
  - `pvdc_stagetwo.model`
- The training and testing datasets used to fit and validate the MLIP model:
  - `training-set.xyz`
  - `testing-set.xyz`

A subset of the training and testing datasets were relabelled from the following open sources:

> **Amorphous carbon and graphitic configurations**: Rowe, P.; Deringer, V. L.; Gasparotto, P.; Csányi, G.; Michaelides, A. *J. Chem. Phys.* **2020**, *153*, 034702.

> **Bulk CH configurations**: Ibragimova, R.; Kuklin, M. S.; Zarrouk, T.; Caro, M. A. *Chem. Mater.* **2025**, *37*, 1094–1110.

> **Hydrocarbon molecules**: Kim, S.; Chen, J.; Cheng, T.; Gindulyte, A.; He, J.; He, S.; Li, Q.; Shoemaker, B. A.;Thiessen, P. A.; Yu, B.; Zaslavsky, L.; Zhang, J.; Bolton, E. E. *Nucleic Acids Res.* **2025**, *53*, D1516–D1525.

> **QM9 molecules**: Wichmann, C.; Ramakrishnan, R.; von Lilienfeld, O. A. *Sci. Data* **2025**, *12*, 9.

### `scripts/`

The analysis scripts, processed data, and plotting scripts are organized by manuscript section:

- `1a/` — representative trajectory snapshots at multiple time points, including an OVITO session file for visualization
- `1b/` — simulation temperature profile
- `1c/` — mass over time
- `1d/` — density over time
- `2a/` — bond count over time
- `2d/` — pre-reactive encounter time
- `2e/` — event-based lifetime and frequency
- `2f/` — reaction type based on regioselectivity
- `3a/` — carbon coordination and cluster over time
- `3d/` — representative trajectory snapshots at multiple time points with six-membered rings highlighted, together with the corresponding ring-identification script

Each section directory typically contains:

- an analysis script (`*_analysis.py`),
- a plotting script (`*_plot.py`), and
- a generated CSV summary (`*.csv`).
