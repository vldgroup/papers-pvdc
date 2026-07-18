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
OUTPUT_CSV = "event_based_lifetime_freq.csv"

cutoff_c_cl_rad = 2.3
cutoff_c_h_rad = 1.4
cutoff_h_cl_rad = 1.6
cutoff_cl_cl_rad = 2.4
cutoff_h_h_rad = 1.1

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

for frame in range(1, frames):
    positions = data[frame].get_scaled_positions()
    count = len(positions)
    atom_counts.append(count)

    if count == atom_counts[frame - 1]:
        ids.append(list(ids[frame - 1]))
        indices.append(list(indices[frame - 1]))
        continue

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

rad_cutoffs = {
    (6, 17): cutoff_c_cl_rad,
    (17, 6): cutoff_c_cl_rad,
    (6, 1): cutoff_c_h_rad,
    (1, 6): cutoff_c_h_rad,
    (17, 1): cutoff_h_cl_rad,
    (1, 17): cutoff_h_cl_rad,
    (17, 17): cutoff_cl_cl_rad,
    (1, 1): cutoff_h_h_rad,
}

secluded_orig_sets = [set() for _ in range(frames)]

for frame_idx in range(frames):
    atoms = data[frame_idx]
    n_local = len(atoms)
    if n_local == 0:
        continue

    orig_of_local = np.full(n_local, -1, dtype=int)
    for orig_id, local_idx in enumerate(indices[frame_idx]):
        if local_idx != -1 and 0 <= local_idx < n_local:
            orig_of_local[local_idx] = orig_id

    try:
        i_idx, j_idx = neighbor_list("ij", atoms, cutoff=rad_cutoffs)
    except Exception:
        i_idx = np.array([], dtype=int)
        j_idx = np.array([], dtype=int)

    has_radial_neighbor = np.zeros(n_local, dtype=bool)
    if i_idx.size:
        has_radial_neighbor[i_idx] = True
        has_radial_neighbor[j_idx] = True

    secluded_local = ~has_radial_neighbor

    for local_idx in np.nonzero(secluded_local)[0]:
        orig_id = int(orig_of_local[local_idx])
        if orig_id < 0:
            continue

        symbol = atoms.get_chemical_symbols()[local_idx]
        if symbol in ("H", "Cl"):
            secluded_orig_sets[frame_idx].add(orig_id)

rows_written = 0

with open(OUTPUT_CSV, "w", newline="") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(
        ["Symbol", "Atom_ID", "Event_index", "Duration_ps", "Start_time_ps", "End_time_ps"]
    )

    for orig_id in range(atom_counts[0]):
        symbol = None
        for frame_idx in range(frames):
            local_idx = indices[frame_idx][orig_id]
            if local_idx != -1:
                symbol = data[frame_idx].get_chemical_symbols()[local_idx]
                break

        if symbol not in ("H", "Cl"):
            continue

        event_count = 0
        frame_idx = 0

        while frame_idx < frames:
            if orig_id not in secluded_orig_sets[frame_idx]:
                frame_idx += 1
                continue

            start_frame = frame_idx
            while frame_idx + 1 < frames and orig_id in secluded_orig_sets[frame_idx + 1]:
                frame_idx += 1
            end_frame = frame_idx

            start_time = frame_times[start_frame - 1] if start_frame > 0 else 0.0
            end_time = frame_times[end_frame]
            duration = end_time - start_time

            if duration > 0:
                event_count += 1
                writer.writerow(
                    [
                        symbol,
                        orig_id + 1,
                        event_count,
                        round(duration, 4),
                        round(start_time, 4),
                        round(end_time, 4),
                    ]
                )
                rows_written += 1

            frame_idx = end_frame + 1

print(f"Saved {rows_written} rows to {OUTPUT_CSV}")
