import os

import numpy as np

from ase import Atoms
from ase.calculators.vasp import Vasp

from ase.io import read


if __name__ == '__main__':


    npar = int(os.environ["SLURM_NNODES"])

    arguments = {"non_magnetic" : {"ispin": 1},
                 "magnetic_001" : {"ispin": 2, # spin polarised calculation
                                   "magmom" : [0.0, -3.0, 3.0]}, # a first guess for the Mn2}
                 "magnetic_1-10" : {"ispin": 2,
                                   "lnoncollinear": True,
                                   "magmom" : [[0.0, 0.0, 0.0],
                                              [-3.0 / np.sqrt(2.0), 3.0 / np.sqrt(2.0), 0.0],
                                               [3.0 / np.sqrt(2.0), -3.0 / np.sqrt(2.0), 0.0]]},
                 "magnetic_110" : {"ispin": 2,
                                   "lnoncollinear": True,
                                   "magmom" : [[0.0, 0.0, 0.0],
                                               [3.0 / np.sqrt(2.0), 3.0 / np.sqrt(2.0), 0.0],
                                               [-3.0 / np.sqrt(2.0), -3.0 / np.sqrt(2.0), 0.0]]}}

    for name, kwargs in arguments.items():
        print(kwargs)
        VaspCalc = Vasp(npar=npar,
                        lorbit=11,
                        kpts=[11, 11, 11],
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
                        xc='PBE',
                        **kwargs)

        primitive_unit_cell = read("strain_filtered_unit_cell.xyz")
        primitive_unit_cell.calc = VaspCalc
        energy = primitive_unit_cell.get_potential_energy()
        print(f"Evaluation results with {name}")
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
        primitive_unit_cell.write(f"reevaluated_primitive_cell_{name}.xyz")
