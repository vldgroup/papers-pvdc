import numpy as np
import pandas as pd
from ase.io import read

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
OUTPUT_CSV = "mass_over_time.csv"

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
    masses_da = np.array([float(np.sum(atoms.get_masses())) for atoms in frames])
    frame_idx = np.arange(len(masses_da))
    time_ps = t_start + frame_idx * dt

    for frame, t_ps, mass_da in zip(frame_idx, time_ps, masses_da):
        rows.append(
            {
                "stage": stage,
                "frame": frame,
                "time_ps": t_ps,
                "total_mass_da": mass_da,
            }
        )

df = pd.DataFrame(rows)
df = df.sort_values("time_ps").reset_index(drop=True)
df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved {len(df)} rows to {OUTPUT_CSV}")
