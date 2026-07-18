import networkx as nx
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
OUTPUT_CSV = "c_cn_cluster_over_time.csv"
CUTOFFS = {(6, 1): 1.30, (6, 6): 1.85, (6, 17): 2.10}
CUTOFF_CC = 1.85
N_CARBON_TOTAL = 1176

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

    largest_cluster_sizes = []
    cn4_counts = []
    cn3_counts = []
    cn2_counts = []
    n_carbon_per_frame = []

    for atoms in frames:
        c_indices = [
            i for i, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == "C"
        ]
        n_carbon_frame = len(c_indices)
        n_carbon_per_frame.append(n_carbon_frame)

        c_atoms = atoms[c_indices]
        i_idx, j_idx = neighbor_list(
            "ij", c_atoms, cutoff=CUTOFF_CC, self_interaction=False
        )

        graph = nx.Graph()
        graph.add_nodes_from(range(n_carbon_frame))
        mask = i_idx < j_idx
        edges = (
            np.column_stack((i_idx[mask], j_idx[mask]))
            if np.any(mask)
            else np.empty((0, 2), dtype=int)
        )
        graph.add_edges_from(edges)

        if graph.number_of_edges() > 0:
            largest_cc = max(nx.connected_components(graph), key=len)
            largest_cluster_sizes.append(len(largest_cc))
        else:
            largest_cluster_sizes.append(1 if n_carbon_frame > 0 else 0)

        i_all, _ = neighbor_list("ij", atoms, cutoff=CUTOFFS)
        coordinations = np.bincount(i_all, minlength=len(atoms))
        carbon_coords = coordinations[c_indices]

        cn4_counts.append(np.sum(carbon_coords == 4))
        cn3_counts.append(np.sum(carbon_coords == 3))
        cn2_counts.append(np.sum(carbon_coords == 2))

    n_carbon_per_frame = np.array(n_carbon_per_frame)
    largest_cluster_sizes = np.array(largest_cluster_sizes)
    cn4_counts = np.array(cn4_counts)
    cn3_counts = np.array(cn3_counts)
    cn2_counts = np.array(cn2_counts)

    for (
        frame,
        t_ps,
        cluster_size,
        cn4_count,
        cn3_count,
        cn2_count,
        n_carbon_frame,
    ) in zip(
        frame_idx,
        time_ps,
        largest_cluster_sizes,
        cn4_counts,
        cn3_counts,
        cn2_counts,
        n_carbon_per_frame,
    ):
        rows.append(
            {
                "stage": stage,
                "frame": frame,
                "time_ps": t_ps,
                "largest_cluster_size": cluster_size,
                "largest_cluster_pct": (cluster_size / N_CARBON_TOTAL) * 100.0,
                "cn4_pct": (
                    (cn4_count / n_carbon_frame) * 100.0 if n_carbon_frame else 0.0
                ),
                "cn3_pct": (
                    (cn3_count / n_carbon_frame) * 100.0 if n_carbon_frame else 0.0
                ),
                "cn2_pct": (
                    (cn2_count / n_carbon_frame) * 100.0 if n_carbon_frame else 0.0
                ),
            }
        )

df = pd.DataFrame(rows)
df = df.sort_values("time_ps").reset_index(drop=True)
df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved {len(df)} rows to {OUTPUT_CSV}")
