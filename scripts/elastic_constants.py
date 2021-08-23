import os
import numpy as np

from ase import Atoms
from ase.io import read
from ase.calculators.vasp import Vasp

from ase.constraints import StrainFilter, UnitCellFilter, ExpCellFilter
from ase.optimize import FIRE, QuasiNewton
from matscipy.elasticity import fit_elastic_constants

def make_primitive_Mn2Au_cell(alat, factor):

    cell = [[-alat / 2.0, alat / 2.0, alat * factor / 2.0],
            [alat / 2.0, -alat / 2.0, alat * factor / 2.0],
            [alat / 2.0, alat / 2.0, -alat * factor / 2.0]
           ]

    Au = Atoms("Au", positions=[[0.0, 0.0, 0.0]], cell=cell, pbc=True)
    Mn2 = Atoms("Mn2", scaled_positions=[[1.0/ 3.0, 1.0 / 3.0, 0.0], [2.0 / 3.0, 2.0 / 3.0, 0.0]], cell=cell, pbc=True)

    return Au + Mn2

if __name__ == '__main__':

    # initial guesses from the paper
    alat = 3.28013797
    factor = 2.5713

    npar = int(os.environ["SLURM_NNODES"])

    VaspCalc = Vasp(npar=npar,
                    kpar=10,
                    lorbit=11,
                    kpts=[25, 25, 25],
                    istart=1, # read wavecar
                    # icharg=1, # read chgcar
                    lwave=False, # do not write wavecar
                    lcharg=False, # do not write CHGCAR
                    # istart=0, # start from scratch
                    # icharg=2, # default for istart=0
                    isif=2,
                    # nsim=2,
                    prec='Accurate',
                    encut=600,
                    ediff=1.e-7,
                    nelm=200,
                    nelmin=5,
                    algo="Normal",
                    isym=-1,
                    ismear=1,
                    sigma=0.1,
                    ispin=2, # spin polarised calculation
                    magmom=[0.0, -3.5, 3.5], # a first guess for the Mn2Au
                    xc='PBE')

    # primitive_unit_cell = make_primitive_Mn2Au_cell(alat, factor)
    primitive_unit_cell = read("strain_filtered_unit_cell.xyz")
    primitive_unit_cell.calc = VaspCalc

    print("Before strain filter")
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
    fmax = np.linalg.norm(primitive_unit_cell.get_forces(), axis=1).max()
    print(f"Fmax = {fmax} eV")
    # sf = StrainFilter(primitive_unit_cell)

    # opt = FIRE(sf)
    # opt.run(fmax=1e-4, steps=10)  # max force in eV/A
    
    # opt = QuasiNewton(primitive_unit_cell)
    # opt.run(fmax=1.0e-3)i
    
    constants, errors = fit_elastic_constants(primitive_unit_cell, symmetry="tetragonal_low", verbose=True)
    
    print(constants)
    print(errors)

    np.savetxt("constant.csv", constants)
    np.savetxt("errors.csv", errors)

    print("After force_min filter")
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
    fmax = np.linalg.norm(primitive_unit_cell.get_forces(), axis=1).max()
    print(f"Fmax = {fmax} eV")

    primitive_unit_cell.write("strain_filtered_unit_cell_forces.xyz")
