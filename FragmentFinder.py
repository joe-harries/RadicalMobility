import numpy as np
from scipy.spatial.distance import pdist, squareform

def read_xyz(file_path):
    atoms = []
    coordinates = []

    with open(file_path, 'r') as file:
        lines = file.readlines()[2:]  # Skip first two lines (atom count and comment)
        for line in lines:
            parts = line.split()
            atom = parts[0]
            x, y, z = map(float, parts[1:])
            atoms.append(atom)
            coordinates.append([x, y, z])
    
    return atoms, np.array(coordinates)

def identify_fragments(atoms, coordinates, bond_threshold=1.5):
    distance_matrix = squareform(pdist(coordinates))
    n_atoms = len(atoms)
    visited = [False] * n_atoms
    fragments = []

    def dfs(atom_index, fragment):
        visited[atom_index] = True
        fragment.append(atom_index)
        for j in range(n_atoms):
            if not visited[j] and distance_matrix[atom_index][j] < bond_threshold:
                dfs(j, fragment)

    for i in range(n_atoms):
        if not visited[i]:
            fragment = []
            dfs(i, fragment)
            fragments.append(fragment)

    return fragments

def print_fragments(fragments, atoms):
    for i, fragment in enumerate(fragments):
        fragment_atoms = [f"{atoms[idx]}({idx + 1})" for idx in fragment]  # +1 for 1-based index
        print(f"Fragment {i + 1}: {', '.join(fragment_atoms)}")

# Usage
xyz_file = 'untitled.xyz'  # Replace with your actual .xyz file path
atoms, coordinates = read_xyz(xyz_file)
fragments = identify_fragments(atoms, coordinates)
print_fragments(fragments, atoms)
