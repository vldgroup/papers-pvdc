import csv
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
OUTPUT_CSV = "reaction_type.csv"

cutoff_c_cl = 2.1
cutoff_c_h = 1.3
cutoff_c_c = 1.85
cutoff_h_cl = 1.3
cutoff_h_cl_loose = 1.5

adj_cache = {}


def get_adj_atoms(frame, index, cutoff):
    # Cache neighbor lists per frame and cutoff to avoid recomputing
    key = (frame, cutoff)
    if key not in adj_cache:
        atoms = data[frame]
        i_idx, j_idx = neighbor_list("ij", atoms, cutoff=cutoff)
        adj = [[] for _ in range(len(atoms))]
        for i, j in zip(i_idx, j_idx):
            adj[i].append(j)
            adj[j].append(i)
        adj_cache[key] = adj
    return adj_cache[key][index]


def get_carbon_distance(frame, index1, index2, max_distance=3):
    visited = set()
    queue = [(index1, 0)]
    symbols = data[frame].get_chemical_symbols()

    while queue:
        current_index, distance = queue.pop(0)
        if current_index == index2:
            return distance
        if distance >= max_distance:
            continue

        visited.add(current_index)
        adj_atoms = get_adj_atoms(frame, current_index, cutoff_c_c)
        for adj_index in adj_atoms:
            if symbols[adj_index] != "C":
                continue
            if adj_index not in visited:
                queue.append((adj_index, distance + 1))

    return -1


def trace_parent_carbon(atom_id, start_frame, cutoff):
    # Walk backwards from start_frame to the most recent frame in which the
    # atom is bonded to a carbon; return that carbon's ID (-1 if none).
    for f in range(start_frame, -1, -1):
        idx_f = indices[f][atom_id]
        if idx_f == -1:
            continue
        symbols_f = data[f].get_chemical_symbols()
        for adj_idx in get_adj_atoms(f, idx_f, cutoff):
            if symbols_f[adj_idx] == "C":
                return ids[f][adj_idx]
    return -1


data = []
frame_times = []
current_total_time = 0.0

for traj_path, dt in zip(TRAJ_FILES, TIMESTEP_PS):
    traj_data = read(
        traj_path,
        index=":",
        format="lammps-dump-text",
        specorder=["C", "Cl", "H"],
    )
    data.extend(traj_data)
    for _ in range(len(traj_data)):
        frame_times.append(current_total_time)
        current_total_time += dt

frames = len(data)
if frames == 0:
    raise RuntimeError("No trajectory frames were loaded.")

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
        continue

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


last_carbon_ids = [-1] * atom_counts[0]
last_carbon_frames = [-1] * atom_counts[0]
first_h_or_cl_ids = [-1] * atom_counts[0]
first_h_or_cl_frames = [-1] * atom_counts[0]
eliminated_cl_ids = []

for frame in key_frames:
    print(f"Frame {frame - 1} -> {frame} ({frame_times[frame]:.3f} ps)")

    for atom_id in range(atom_counts[0]):
        if indices[frame][atom_id] == -1 and indices[frame - 1][atom_id] != -1:
            index_prev = indices[frame - 1][atom_id]
            symbol = data[frame - 1].get_chemical_symbols()[index_prev]

            if symbol == "Cl":
                eliminated_cl_ids.append(atom_id)

            complement_symbol = (
                "H" if symbol == "Cl" else "Cl" if symbol == "H" else None
            )

            if complement_symbol is None:
                continue

            cutoff_c_x = cutoff_c_cl if symbol == "Cl" else cutoff_c_h
            last_h_or_cl_id = -1
            adj_atoms_elim = get_adj_atoms(frame - 1, index_prev, cutoff_h_cl)
            for adj_index in adj_atoms_elim:
                if (
                    data[frame - 1].get_chemical_symbols()[adj_index]
                    == complement_symbol
                ):
                    last_h_or_cl_id = ids[frame - 1][adj_index]
                    break

            found_partner_frame = frame - 1
            first_encounter_partner_id = last_h_or_cl_id
            last_bonded_carbon_frame = 0
            found_last_c = -1

            for f in range(frame - 1, -1, -1):
                idx_f = indices[f][atom_id]
                if idx_f == -1:
                    continue

                nbrs_bonded = get_adj_atoms(f, idx_f, cutoff_c_x)
                has_carbon_bond = False
                for adj_idx in nbrs_bonded:
                    if data[f].get_chemical_symbols()[adj_idx] == "C":
                        last_bonded_carbon_frame = f
                        found_last_c = ids[f][adj_idx]
                        has_carbon_bond = True
                        break

                adj_loose = get_adj_atoms(f, idx_f, cutoff_h_cl_loose)
                for adj_idx in adj_loose:
                    if data[f].get_chemical_symbols()[adj_idx] == complement_symbol:
                        found_partner_frame = f
                        first_encounter_partner_id = ids[f][adj_idx]
                        break

                if has_carbon_bond:
                    break

            last_carbon_ids[atom_id] = found_last_c
            last_carbon_frames[atom_id] = last_bonded_carbon_frame
            first_h_or_cl_ids[atom_id] = first_encounter_partner_id
            first_h_or_cl_frames[atom_id] = found_partner_frame

print(
    f"Tracking phase complete: {len(key_frames)} key frames, "
    f"{len(eliminated_cl_ids)} eliminated Cl candidates."
)
print(f"Starting CSV write phase -> {OUTPUT_CSV}")

csv_struct_file = open(OUTPUT_CSV, "w", newline="")
csv_struct_writer = csv.writer(csv_struct_file)
csv_struct_writer.writerow(
    [
        "Cl_ID",
        "H_ID",
        "Connectivity",
        "C_Cl_index",
        "C_H_index",
        "C_C_euclidean_distance",
        "Reaction_frame",
        "Time_ps",
        "H_parent_traced_from",
    ]
)

for atom_id in eliminated_cl_ids:
    if first_h_or_cl_ids[atom_id] == -1:
        continue

    frame = last_carbon_frames[atom_id]
    carbon_cl_id = last_carbon_ids[atom_id]

    partner_h_id = first_h_or_cl_ids[atom_id]
    if partner_h_id < 0 or partner_h_id >= len(last_carbon_ids):
        continue
    carbon_h_id = last_carbon_ids[partner_h_id]
    h_parent_source = "elimination"

    # The reacting H was never eliminated, so it was not traced from a removal
    # frame. Trace it back from the frame in which the Cl first met it instead.
    if carbon_h_id == -1:
        carbon_h_id = trace_parent_carbon(
            partner_h_id, first_h_or_cl_frames[atom_id], cutoff_c_h
        )
        h_parent_source = "encounter"

    if frame == -1 or carbon_cl_id == -1 or carbon_h_id == -1:
        continue

    carbon_cl_index = indices[frame][carbon_cl_id]
    carbon_h_index = indices[frame][carbon_h_id]

    carbon_distance = get_carbon_distance(
        frame, carbon_cl_index, carbon_h_index, max_distance=3
    )

    if carbon_distance == 0:
        distance_str = "1,1"
    elif carbon_distance == 1:
        distance_str = "1,2"
    elif carbon_distance == 2:
        distance_str = "1,3"
    elif carbon_distance == 3:
        distance_str = "1,4"
    else:
        distance_str = "others"

    # Use ASE's get_distance with mic=True to handle PBC minimum-image convention
    c_c_euclidean_dist = round(
        float(data[frame].get_distance(carbon_cl_index, carbon_h_index, mic=True)), 4
    )

    csv_struct_writer.writerow(
        [
            atom_id + 1,
            partner_h_id + 1,
            distance_str,
            carbon_cl_index + 1,
            carbon_h_index + 1,
            c_c_euclidean_dist,
            frame,
            round(frame_times[frame], 4),
            h_parent_source,
        ]
    )

csv_struct_file.close()
print(f"Wrote {OUTPUT_CSV}")
