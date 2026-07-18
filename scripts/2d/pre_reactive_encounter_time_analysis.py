import pandas as pd
import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list

TRAJ_FILES = [
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/1_3528_npt_ramp_1200.pos",
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/2_3528_nvt_ramp_3000.pos",
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/3_3528_nvt_hold_3000.pos",
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/4_3528_npt_quench_300.pos",
]

TIMESTEP_PS = [0.1, 0.01, 0.01, 0.1]
OUTPUT_CSV = "pre_reactive_encounter_time.csv"


all_data = []
frame_times = []
current_total_time = 0.0

for traj_path, dt in zip(TRAJ_FILES, TIMESTEP_PS):
    traj_data = read(
        traj_path, index=":", format="lammps-dump-text", specorder=["C", "Cl", "H"]
    )
    all_data.extend(traj_data)

    for _ in range(len(traj_data)):
        frame_times.append(current_total_time)
        current_total_time += dt

data = all_data
frames = len(data)

atom_counts = [len(data[0].get_scaled_positions())]
ids = [list(range(atom_counts[0]))]
indices = [list(range(atom_counts[0]))]
key_frames = []

for frame in range(1, frames):
    positions = data[frame].get_scaled_positions()
    count = len(positions)
    atom_counts.append(count)

    if count == atom_counts[frame - 1]:
        ids.append(list(ids[frame - 1]))
        indices.append(list(indices[frame - 1]))
    else:
        key_frames.append(frame)
        prev_positions = data[frame - 1].get_scaled_positions()
        new_ids = []

        for i in range(count):
            matched_id = -1
            for prev_idx, prev_pos in enumerate(prev_positions):
                delta = positions[i] - prev_pos
                delta -= np.rint(delta)
                if np.linalg.norm(delta) < 1e-3:
                    matched_id = ids[frame - 1][prev_idx]
                    break
            new_ids.append(matched_id)

        ids.append(new_ids)

        new_indices = [-1] * atom_counts[0]
        for idx, atom_id in enumerate(new_ids):
            if atom_id != -1:
                new_indices[atom_id] = idx
        indices.append(new_indices)

cutoff_c_cl, cutoff_c_h, cutoff_h_cl = 2.1, 1.3, 1.3
cutoff_h_cl_loose = 1.5

cached_neighbors_loose = {}
cached_neighbors_hcl = {}
cached_neighbors_c_cl = {}
cached_neighbors_c_h = {}

for frame_idx in range(frames):
    atoms = data[frame_idx]

    try:
        i_idx, j_idx = neighbor_list("ij", atoms, cutoff=cutoff_h_cl_loose)
    except Exception:
        i_idx = j_idx = []
    neighbors_loose = {}
    for i, j in zip(i_idx, j_idx):
        neighbors_loose.setdefault(i, []).append(j)
        neighbors_loose.setdefault(j, []).append(i)
    cached_neighbors_loose[frame_idx] = neighbors_loose

    try:
        i_idx, j_idx = neighbor_list("ij", atoms, cutoff=cutoff_h_cl)
    except Exception:
        i_idx = j_idx = []
    neighbors_hcl = {}
    for i, j in zip(i_idx, j_idx):
        neighbors_hcl.setdefault(i, []).append(j)
        neighbors_hcl.setdefault(j, []).append(i)
    cached_neighbors_hcl[frame_idx] = neighbors_hcl

    try:
        i_idx, j_idx = neighbor_list("ij", atoms, cutoff=cutoff_c_cl)
    except Exception:
        i_idx = j_idx = []
    neighbors_c_cl = {}
    for i, j in zip(i_idx, j_idx):
        neighbors_c_cl.setdefault(i, []).append(j)
        neighbors_c_cl.setdefault(j, []).append(i)
    cached_neighbors_c_cl[frame_idx] = neighbors_c_cl

    try:
        i_idx, j_idx = neighbor_list("ij", atoms, cutoff=cutoff_c_h)
    except Exception:
        i_idx = j_idx = []
    neighbors_c_h = {}
    for i, j in zip(i_idx, j_idx):
        neighbors_c_h.setdefault(i, []).append(j)
        neighbors_c_h.setdefault(j, []).append(i)
    cached_neighbors_c_h[frame_idx] = neighbors_c_h

rows = []

for frame in key_frames:
    for atom_id in range(atom_counts[0]):
        if indices[frame][atom_id] == -1 and indices[frame - 1][atom_id] != -1:
            index_prev = indices[frame - 1][atom_id]
            symbol = data[frame - 1].get_chemical_symbols()[index_prev]

            complement_symbol = (
                "H" if symbol == "Cl" else "Cl" if symbol == "H" else None
            )
            if complement_symbol is None:
                continue

            found_partner_frame = frame - 1
            last_bonded_carbon_frame = 0
            first_encounter_partner_id = -1

            for f in range(frame - 1, -1, -1):
                idx_f = indices[f][atom_id]
                if idx_f == -1:
                    continue

                atoms_f = data[f]
                symbols_f = atoms_f.get_chemical_symbols()

                has_carbon_bond = False
                if symbol == "Cl":
                    neighbors_c = cached_neighbors_c_cl.get(f, {})
                else:
                    neighbors_c = cached_neighbors_c_h.get(f, {})

                if idx_f in neighbors_c:
                    for adj_idx in neighbors_c[idx_f]:
                        if symbols_f[adj_idx] == "C":
                            last_bonded_carbon_frame = f
                            has_carbon_bond = True
                            break

                neighbors_loose_f = cached_neighbors_loose.get(f, {})
                if not has_carbon_bond and idx_f in neighbors_loose_f:
                    for adj_idx in neighbors_loose_f[idx_f]:
                        if symbols_f[adj_idx] == complement_symbol:
                            found_partner_frame = f
                            first_encounter_partner_id = ids[f][adj_idx]
                            break

                if has_carbon_bond:
                    break

            index_prev = indices[frame - 1][atom_id]
            partner_id_before = -1
            if index_prev != -1:
                neighbors_prev = cached_neighbors_hcl.get(frame - 1, {})
                if index_prev in neighbors_prev:
                    for adj_idx in neighbors_prev[index_prev]:
                        if (
                            data[frame - 1].get_chemical_symbols()[adj_idx]
                            == complement_symbol
                        ):
                            partner_id_before = ids[frame - 1][adj_idx]
                            break

            t_1st_iso = (
                frame_times[found_partner_frame] - frame_times[last_bonded_carbon_frame]
            )

            atom_info = f"{atom_id + 1} ({index_prev + 1} @ {frame - 1})"
            rows.append(
                {
                    "Symbol": symbol,
                    "Atom_ID": atom_info,
                    "Partner_ID": (
                        partner_id_before + 1 if partner_id_before != -1 else -1
                    ),
                    "First_partner_ID": (
                        first_encounter_partner_id + 1
                        if first_encounter_partner_id != -1
                        else -1
                    ),
                    "1st_iso_time_ps": round(t_1st_iso, 4),
                }
            )


df = pd.DataFrame(rows)
df = df[df["Symbol"].isin(["H", "Cl"])].reset_index(drop=True)
df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved {len(df)} rows to {OUTPUT_CSV}")
