import os
from optimize_geometry import optimize_for_states, optimize_geometry
from reorg_energy import calculate_reorganization_energy
from stack_dimer import stack_dimer
from TransferIntegral import calculate_transfer_integral
from mobility import calculate_charge_carrier_mobility

log_files = {
    'radical_opt': 'radical_optimization.log',
    'anion_opt': 'anion_optimization.log',
    'cation_opt': 'cation_optimization.log',
    'reorg_energy': 'reorganization_energy.log',
    'stack_opt': 'stack_optimization.log',
    'transfer_integral': 'transfer_integral.log',
    'mobility_calc': 'mobility_calculation.log'
}

def clear_restart_files():
    restart_files = ['wbo', 'xtbrestart', 'xtbopt.mol', 'xtbopt.xyz']
    for file in restart_files:
        if os.path.exists(file):
            os.remove(file)

def main():
    optimized_files = optimize_for_states('Blatter1.mol')
    optimized_radical_file = optimized_files['radical']
    optimized_anion_file = optimized_files['anion']
    optimized_cation_file = optimized_files['cation']

    reorg_energy_electron, reorg_energy_hole = calculate_reorganization_energy(
        optimized_radical_file, optimized_anion_file, optimized_cation_file
    )

    dimer_file = stack_dimer(optimized_radical_file, 3.5)  
    optimized_dimer_file = "optimized_dimer.xyz"
    clear_restart_files()
    optimize_geometry(dimer_file, optimized_dimer_file, 0, 2)

    transfer_integral_electron = calculate_transfer_integral(optimized_dimer_file)
    transfer_integral_hole = calculate_transfer_integral(optimized_dimer_file)

    electron_mobility = calculate_charge_carrier_mobility(
        reorg_energy_electron, transfer_integral_electron, temperature=300
    )
    hole_mobility = calculate_charge_carrier_mobility(
        reorg_energy_hole, transfer_integral_hole, temperature=300
    )

    print("\n--- Calculation Results ---")
    print(f"Reorganization Energy (Electron): {reorg_energy_electron:.4f} eV")
    print(f"Reorganization Energy (Hole): {reorg_energy_hole:.4f} eV")
    print(f"Transfer Integral (Electron): {transfer_integral_electron:.4f} eV")
    print(f"Transfer Integral (Hole): {transfer_integral_hole:.4f} eV")
    print(f"Electron Mobility: {electron_mobility:.4f} cm²/V·s")
    print(f"Hole Mobility: {hole_mobility:.4f} cm²/V·s")

if __name__ == "__main__":
    main()
