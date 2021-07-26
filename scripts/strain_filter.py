import os

from ase import Atoms
from ase.calculators.vasp import Vasp

from ase.constraints import StrainFilter

def make_primitive_Mn2Au_cell(alat, factor):

    cell = [[-alat / 2.0, alat / 2.0, alat * factor / 2.0],
            [alat / 2.0, -alat / 2.0, alat * factor / 2.0],
            [alat / 2.0, alat / 2.0, -alat * factor / 2.0]
           ]

    Au = Atoms("Au", positions=[[0.0, 0.0, 0.0]], cell=cell, pbc=True)
    Mn2 = Atoms("Mn2", scaled_positions=[[0.333, 0.333, 0.0], [0.6666, 0.666, 0.0]], cell=cell, pbc=True)

    return Au + Mn2

if __name__ == '__main__':

    # initial guesses from the paper    
    alat = 3.1445
    factor = 2.5658

    npar = int(os.environ["SLURM_NNODES"])

    VaspCalc = Vasp(npar=npar,
                    lorbit=11,
                    kpts=[6, 6, 6],
                    istart=0, # start from scratch
                    icharg=2, # default for istart=0
                    isif=2,
                    nsim=2,
                    prec='Accurate',
                    encut=500,
                    ediff=1.e-6,
                    nelm=200,
                    nelmin=5,
                    algo="Normal",
                    isym=0,
                    ismear=1,
                    sigma=0.1,
                    ispin=2, # spin polarised calculation
                    magmom=[0.0, -3.0, 3.0], # a first guess for the Mn2Au
                    xc='PBE')

    primitive_unit_cell = make_primitive_Mn2Au_cell(alat, factor)
    primitive_unit_cell.calc = VaspCalc

    print("Before strain filter")
    print(primitive_unit_cell.cell.get_bravais_lattice())

    sf = StrainFilter(primitive_unit_cell)

    opt = FIRE(sf)
    opt.run(fmax=1e-4)  # max force in eV/A

    print("After strain filter")
    print(primitive_unit_cell.cell.get_bravais_lattice())