import os

from ase.lattice.tetragonal import CenteredTetragonalFactory

from ase import Atoms
from ase.calculators.vasp import Vasp

from ase.constraints import StrainFilter, UnitCellFilter, ExpCellFilter
from ase.optimize import FIRE

class Mn2AuFactory(CenteredTetragonalFactory):
    "A factory for creating Mn2Au lattices."
    bravais_basis = [[0.0, 0.0,  0.333333], [0.0 , 0.0, 0.666666], [0, 0, 0]]
    element_basis = (0, 0, 1)

if __name__ == '__main__':

    # initial guesses from the paper
    alat = 3.27923
    factor = 2.5713

    npar = int(os.environ["SLURM_NNODES"])


    Mn2Au = Mn2AuFactory()
    cubic_unit_cell = Mn2Au(["Mn", "Au"], latticeconstant=[alat, factor * alat])
    cubic_unit_cell.wrap()

    VaspCalc = Vasp(npar=npar, kpts=[9, 9, 6], # nupdown=0.0,
                    # amix=0.2,
                    # bmix=0.0001,
                    # amix_mag=0.8,
                    # bmix_mag=0.0001,
                    lorbit=11,
                    istart=0, # start from scratch
                    icharg=2, # default for istart=0
                    isif=2,
                    nsim=2,
                    prec='Accurate',
                    encut=500,
                    ediff=1.e-6,
                    ediffg=-0.04,
                    nelm=100,
                    nelmin=5,
                    algo="Normal",
                    isym=0,
                    ismear=1,
                    sigma=0.1,
                    ispin=2, # spin polarised calculation
                    magmom=[-3.0, 3.0, 3.0, -3.0, 0.0, 0.0], # a first guess for the Mn2Au the order of atoms in VASP is bizarre
                    xc='PBE')

    cubic_unit_cell.calc = VaspCalc
    energy = cubic_unit_cell.get_potential_energy()

    print("Before strain filter")
    print(cubic_unit_cell.cell.get_bravais_lattice())

    energy = cubic_unit_cell.get_potential_energy()
    print(energy)
    print(f"Energy per atom: {energy/len(cubic_unit_cell)}")
    print("Forces")
    print(cubic_unit_cell.get_forces())
    print("Stress")
    print(cubic_unit_cell.get_stress())
    print("Total magnetic moment")
    print(cubic_unit_cell.get_magnetic_moment())
    print("Magnetic moments:")
    print(cubic_unit_cell.get_magnetic_moments())
    print(cubic_unit_cell.cell)
    sf = ExpCellFilter(cubic_unit_cell)

    opt = FIRE(sf)
    opt.run(fmax=1e-3)  # max force in eV/A

    print("After strain filter")
    print(cubic_unit_cell.cell.get_bravais_lattice())

    energy = cubic_unit_cell.get_potential_energy()
    print(energy)
    print(f"Energy per atom: {energy/len(cubic_unit_cell)}")
    print("Forces")
    print(cubic_unit_cell.get_forces())
    print("Stress")
    print(cubic_unit_cell.get_stress())
    print("Total magnetic moment")
    print(cubic_unit_cell.get_magnetic_moment())
    print("Magnetic moments:")
    print(cubic_unit_cell.get_magnetic_moments())

    print(cubic_unit_cell.cell)

    cubic_unit_cell.write("strain_filtered_unit_cell.xyz")
