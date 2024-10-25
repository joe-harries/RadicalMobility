import subprocess

def calculate_single_molecule_energy(optimized_xyz_file, charge, multiplicity):
    result = subprocess.run(['xtb', optimized_xyz_file, '--gfn2', '--sp', '--chrg', str(charge), '--uhf', str(multiplicity - 1)], 
                            check=True, capture_output=True, text=True)
    energy_line = next(line for line in result.stdout.splitlines() if 'TOTAL ENERGY' in line)
    return float(energy_line.split()[3])

def calculate_reorganization_energy(neutral_xyz_file, anion_xyz_file, cation_xyz_file):
    """
    Calculates both electron and hole reorganization energies based on the provided
    neutral, anion, and cation geometries.
    """
    conversion_factor = 27.2114  # Conversion factor from Hartree to eV

    # --- Electron Reorganization Energy Calculation ---
    # Neutral and anion geometries
    E_neutral_optimized = calculate_single_molecule_energy(neutral_xyz_file, 0, 2)  # Neutral geometry optimized
    E_anion_optimized = calculate_single_molecule_energy(anion_xyz_file, -1, 2)  # Anion geometry optimized
    E_neutral_in_anion_geom = calculate_single_molecule_energy(anion_xyz_file, 0, 2)  # Neutral in anion geometry
    E_anion_in_neutral_geom = calculate_single_molecule_energy(neutral_xyz_file, -1, 2)  # Anion in neutral geometry

    # Reorganization energy for electron
    reorg_energy_electron_hartree = (E_neutral_in_anion_geom - E_neutral_optimized) + (E_anion_in_neutral_geom - E_anion_optimized)
    reorg_energy_electron_eV = reorg_energy_electron_hartree * conversion_factor

    # --- Hole Reorganization Energy Calculation ---
    # Neutral and cation geometries
    E_cation_optimized = calculate_single_molecule_energy(cation_xyz_file, 1, 2)  # Cation geometry optimized
    E_neutral_in_cation_geom = calculate_single_molecule_energy(cation_xyz_file, 0, 2)  # Neutral in cation geometry
    E_cation_in_neutral_geom = calculate_single_molecule_energy(neutral_xyz_file, 1, 2)  # Cation in neutral geometry

    # Reorganization energy for hole
    reorg_energy_hole_hartree = (E_neutral_in_cation_geom - E_neutral_optimized) + (E_cation_in_neutral_geom - E_cation_optimized)
    reorg_energy_hole_eV = reorg_energy_hole_hartree * conversion_factor

    return reorg_energy_electron_eV, reorg_energy_hole_eV

if __name__ == "__main__":
    # Input files for the optimized neutral, anion, and cation geometries
    neutral_xyz_file = "optBlatter2.mol"       # Optimized neutral geometry
    anion_xyz_file = "optBlatter2anion.mol"    # Optimized anion geometry
    cation_xyz_file = "optBlatter2cation.mol"  # Optimized cation geometry

    # Call the function to calculate reorganization energies
    reorg_energy_electron, reorg_energy_hole = calculate_reorganization_energy(neutral_xyz_file, anion_xyz_file, cation_xyz_file)

    # Print the results
    print(f"Electron Reorganization Energy: {reorg_energy_electron:.6f} eV")
    print(f"Hole Reorganization Energy: {reorg_energy_hole:.6f} eV")
