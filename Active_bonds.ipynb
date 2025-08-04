# ======== Import Dependencies ========
import numpy as np
import csv
from itertools import combinations

# ======= Read .xyz files ========
def read_xyz(filename):
    """Parse .xyz file and return a list of atoms and a coordinate matrix."""
    with open(filename, 'r') as f:
        lines = f.readlines()
        atom_count = int(lines[0].strip())
        atoms = []
        coordinates = []
        for line in lines[2:2 + atom_count]:
            parts = line.split()
            atoms.append(parts[0])
            coordinates.append([float(x) for x in parts[1:]])
        return atoms, np.array(coordinates)

# ======= Calculate Interatomic Distances =======
def calculate_distance(coord1, coord2):
    """Calculate the Euclidean distance between two 3D points."""
    return np.linalg.norm(coord1 - coord2)

# ======= Identify the list of bonded atom index pairs =======
def find_bonds(atoms, coords, tolerance=0.4):
    """
    Identify bonds based on interatomic distances and covalent radii.
    """
    covalent_radii = {
        "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "F": 0.57
    }
    bonds = []
    atom_count = len(atoms)
    for i, j in combinations(range(atom_count), 2):
        if atoms[i] in covalent_radii and atoms[j] in covalent_radii:
            max_bond_length = covalent_radii[atoms[i]] + covalent_radii[atoms[j]] + tolerance
            distance = calculate_distance(coords[i], coords[j])
            if distance <= max_bond_length:
                bonds.append((i, j))
    return bonds


# ======= Bond lengths comparison across structures (reactant, TS, Product) =======
def compare_bond_lengths(bonds, coords_rea, coords_ts, coords_pro):
    """
    Compare bond lengths for specified bonds across reactant, TS, and product structures.
    """
    results = []
    for atom1, atom2 in bonds:
        bond_rea = calculate_distance(coords_rea[atom1], coords_rea[atom2])
        bond_ts = calculate_distance(coords_ts[atom1], coords_ts[atom2])
        bond_pro = calculate_distance(coords_pro[atom1], coords_pro[atom2])

        delta_rea_ts = bond_ts - bond_rea
        delta_ts_pro = bond_pro - bond_ts

        results.append({
            "atom_pair": (atom1 + 1, atom2 + 1),  # Adjusting to 1-based indexing for readability
            "rea_bond": round(bond_rea, 3),
            "ts_bond": round(bond_ts, 3),
            "pro_bond": round(bond_pro, 3),
            "delta_rea_ts": round(delta_rea_ts, 3),
            "delta_ts_pro": round(delta_ts_pro, 3)
        })
    return results

def filter_significant_changes(results, threshold=0.03):
    """
    Filter bond length changes that are greater than or equal to the given threshold.
    """
    filtered_results = [
        res for res in results
        if abs(res["delta_rea_ts"]) >= threshold or abs(res["delta_ts_pro"]) >= threshold
    ]
    return filtered_results

def save_results_to_csv(results, output_file):
    """Save bond length comparison results to a CSV file."""
    with open(output_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Atom Pair", "REA Bond (Å)", "TS Bond (Å)", "PRO Bond (Å)",
                         "Delta (REA -> TS)", "Delta (TS -> PRO)"])
        for res in results:
            atom1, atom2 = res["atom_pair"]
            writer.writerow([f"{atom1}-{atom2}", res["rea_bond"], res["ts_bond"], res["pro_bond"],
                             res["delta_rea_ts"], res["delta_ts_pro"]])

def main():
    # File Names
    rea_file = "rea.xyz"
    ts_file = "ts.xyz"
    pro_file = "pro.xyz"
    output_file_all = "bond_length_changes.csv"
    output_file_filtered = "significant_bond_changes.csv"

    # Read structures
    atoms_rea, coords_rea = read_xyz(rea_file)
    atoms_ts, coords_ts = read_xyz(ts_file)
    atoms_pro, coords_pro = read_xyz(pro_file)

    # Identify bonds in the TS structure
    bonds = find_bonds(atoms_ts, coords_ts, tolerance=0.4)

    # Compare bond lengths for the identified bonds
    results = compare_bond_lengths(bonds, coords_rea, coords_ts, coords_pro)

    # Save all bond length changes to CSV
    save_results_to_csv(results, output_file_all)
    print(f"All bond length changes saved to {output_file_all}")

    # Filter and save significant bond length changes to CSV
    significant_results = filter_significant_changes(results, threshold=0.03)
    save_results_to_csv(significant_results, output_file_filtered)
    print(f"Significant bond length changes saved to {output_file_filtered}")

if __name__ == "__main__":
    main()
