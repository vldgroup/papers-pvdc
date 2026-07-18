import numpy as np
import pandas as pd
from ase.io import read
from ase.neighborlist import neighbor_list

TRAJ_FILES = [
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/1_3528_npt_ramp_1200.pos",
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/2_3528_nvt_ramp_3000.pos",
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/3_3528_nvt_hold_3000.pos",
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/4_3528_npt_quench_300.pos",
]

STAGE_NAMES = [
    "npt_ramp_1200",
    "nvt_ramp_3000",
    "nvt_hold_3000",
    "npt_quench_300",
]

TIMESTEP_PS = [0.1, 0.01, 0.01, 0.1]
PHASE_STARTS_PS = [0.0, 90.10, 271.90, 372.90]
OUTPUT_CSV = "bonds_counts_over_time.csv"
CUTOFFS = {(6, 1): 1.30, (6, 6): 1.85, (6, 17): 2.10}

rows = []
for stage, traj_path, dt, t_start in zip(
    STAGE_NAMES, TRAJ_FILES, TIMESTEP_PS, PHASE_STARTS_PS
):
    frames = read(
        traj_path,
        index=":",
        format="lammps-dump-text",
        specorder=["C", "Cl", "H"],
    )
    frame_idx = np.arange(len(frames))
    time_ps = t_start + frame_idx * dt

    cc_counts = []
    ch_counts = []
    ccl_counts = []

    for atoms in frames:
        atomic_numbers = atoms.get_atomic_numbers()
        i_idx, j_idx = neighbor_list("ij", atoms, cutoff=CUTOFFS)

        c_mask = atomic_numbers[i_idx] == 6
        j_c = j_idx[c_mask]
        neighbor_z = atomic_numbers[j_c]

        cc_counts.append(np.sum(neighbor_z == 6) / 2)
        ch_counts.append(np.sum(neighbor_z == 1))
        ccl_counts.append(np.sum(neighbor_z == 17))

    for frame, t_ps, cc_count, ch_count, ccl_count in zip(
        frame_idx, time_ps, cc_counts, ch_counts, ccl_counts
    ):
        rows.append(
            {
                "stage": stage,
                "frame": frame,
                "time_ps": t_ps,
                "cc_bonds": cc_count,
                "ch_bonds": ch_count,
                "ccl_bonds": ccl_count,
            }
        )

df = pd.DataFrame(rows)
df = df.sort_values("time_ps").reset_index(drop=True)
df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved {len(df)} rows to {OUTPUT_CSV}")
