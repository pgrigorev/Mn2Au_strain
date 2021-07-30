import os

from ase import Atoms
from ase.calculators.vasp import Vasp


if __name__ == '__main__':

    primitive_unit_cell.read("strain_filtered_unit_cell.xyz")
    print("Old Values")
    print(primitive_unit_cell.cell.get_bravais_lattice())

    energy = primitive_unit_cell.get_potential_energy()
    print(energy)
    print(f"Energy per atom: {energy/len(primitive_unit_cell)}")
    print("Forces")
    print(primitive_unit_cell.get_forces())
    print("Stress")
    print(primitive_unit_cell.get_stress())
    print("Total magnetic moment")
    print(primitive_unit_cell.get_magnetic_moment())
    print("Magnetic moments:")
    print(primitive_unit_cell.get_magnetic_moments())
    print(primitive_unit_cell.cell)

    npar = int(os.environ["SLURM_NNODES"])

    VaspCalc = Vasp(npar=npar,
                    lorbit=11,
                    kpts=[11, 11, 11],
                    istart=1, # read wavecar
                    icharg=2, # default for istart=0
                    isif=2,
                    # nsim=2,
                    prec='Accurate',
                    encut=600,
                    ediff=1.e-6,
                    nelm=200,
                    nelmin=5,
                    algo="Normal",
                    isym=0,
                    ismear=1,
                    sigma=0.1,
                    ispin=2, # spin polarised calculation
                    lsorbit=True, # spin orbit coupling
                    xc='PBE')

    primitive_unit_cell.calc = VaspCalc

    print("Spin-Orbit values:")
    print(primitive_unit_cell.cell.get_bravais_lattice())

    energy = primitive_unit_cell.get_potential_energy()
    print(energy)
    print(f"Energy per atom: {energy/len(primitive_unit_cell)}")
    print("Forces")
    print(primitive_unit_cell.get_forces())
    print("Stress")
    print(primitive_unit_cell.get_stress())
    print("Total magnetic moment")
    print(primitive_unit_cell.get_magnetic_moment())
    print("Magnetic moments:")
    print(primitive_unit_cell.get_magnetic_moments())

    print(primitive_unit_cell.cell)

    primitive_unit_cell.write("strain_filtered_unit_cell.xyz")
