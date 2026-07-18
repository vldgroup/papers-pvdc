import networkx as nx
import numpy as np
from ase.io import read, write
from ase.neighborlist import neighbor_list

INPUT_XYZ = "path/to/your/input_file.xyz"
CUTOFF_CC = 1.85
RING_SIZE = 6
OUTPUT_FILE = "path/to/your/output_file.xyz"

atoms = read(INPUT_XYZ)

ring_flag = np.zeros(len(atoms), dtype=np.int8)
symbols = atoms.get_chemical_symbols()
c_indices = [i for i, symbol in enumerate(symbols) if symbol == "C"]

if len(c_indices) >= RING_SIZE:
    i_idx, j_idx = neighbor_list("ij", atoms, cutoff=CUTOFF_CC, self_interaction=False)

    if i_idx.size:
        c_mask = np.isin(i_idx, c_indices) & np.isin(j_idx, c_indices)
        if np.any(c_mask):
            i_cc = i_idx[c_mask]
            j_cc = j_idx[c_mask]

            global_to_local = {
                global_index: local_index
                for local_index, global_index in enumerate(c_indices)
            }
            edges_local = []

            for a, b in zip(i_cc, j_cc):
                local_a = global_to_local[int(a)]
                local_b = global_to_local[int(b)]
                if local_a < local_b:
                    edges_local.append((local_a, local_b))

            if edges_local:
                graph = nx.Graph()
                graph.add_nodes_from(range(len(c_indices)))
                graph.add_edges_from(edges_local)

                for component in nx.connected_components(graph):
                    if len(component) < RING_SIZE:
                        continue

                    subgraph = graph.subgraph(component)
                    for cycle in nx.minimum_cycle_basis(subgraph):
                        if len(cycle) == RING_SIZE:
                            for local_index in cycle:
                                ring_flag[c_indices[local_index]] = 1

atoms.set_array("ring", ring_flag)
write(OUTPUT_FILE, atoms, format="extxyz")
print(f"Saved {OUTPUT_FILE}")
