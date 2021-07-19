import os

from ase import Atoms
from ase.calculators.vasp import Vasp

from ase.lattice.cubic import SimpleCubicFactory, BodyCenteredCubicFactory

from ase.lattice.tetragonal import CenteredTetragonalFactory

class Mn2AuFactory(CenteredTetragonalFactory):
    "A factory for creating Mn2Au lattices."
    bravais_basis = [[0.0, 0.0,  0.333333], [0.0 , 0.0, 0.666666], [0, 0, 0]]
    element_basis = (0, 0, 1)

if __name__ == '__main__':

    alat = 3.1445
    factor = 2.5658

    npar = int(os.environ["SLURM_NNODES"])


    Mn2Au = Mn2AuFactory()
    cubic_unit_cell = Mn2Au(["Mn", "Au"], latticeconstant=[alat, factor * alat])
    cubic_unit_cell.wrap()


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
                    magmom=[0.0, 2.0, -2.0, 0.0, 2.0, -2.0], # a first guess for the Mn2Au
                    xc='PBE')

    cubic_unit_cell.calc = VaspCalc
    energy = cubic_unit_cell.get_potential_energy()
    print(f"Cubic cell enrgy: {energy}")
