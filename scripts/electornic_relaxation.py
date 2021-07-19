import os

from ase import Atoms
from ase.calculators.vasp import Vasp

def make_primitive_Mn2Au_cell(alat, factor):

    cell = [[-alat / 2.0, alat / 2.0, alat * factor / 2.0],
            [alat / 2.0, -alat / 2.0, alat * factor / 2.0],
            [alat / 2.0, alat / 2.0, -alat * factor / 2.0]
           ]

    Au = Atoms("Au", positions=[[0.0, 0.0, 0.0]], cell=cell, pbc=True)
    Mn2 = Atoms("Mn2", positions=[[0.333, 0.333, 0.0], [0.6666, 0.666, 0.0]], cell=cell, pbc=True)

    return Au + Mn2

if __name__ == '__main__':

    alat = 3.1445
    factor = 2.5658

    npar = int(os.environ["SLURM_NNODES"])

    VaspCalc = Vasp(npar=n_par, kpts=[3,3,3],
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
                    magmom=[0.0, 2.0, -2.0], # a first guess for the Mn2Au
                    xc='PBE')

    primitive_unit_cell = make_primitive_Mn2Au_cell(alat, factor)
    primitive_unit_cell.calc = VaspCalc
    energy = primitive_unit_cell.get_potential_energy()
    print(energy)
