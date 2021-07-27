import os

import pandas as pd

from ase import Atoms
from ase.calculators.vasp import Vasp

def make_primitive_Mn2Au_cell(alat, factor):

    cell = [[-alat / 2.0, alat / 2.0, alat * factor / 2.0],
            [alat / 2.0, -alat / 2.0, alat * factor / 2.0],
            [alat / 2.0, alat / 2.0, -alat * factor / 2.0]
           ]

    Au = Atoms("Au", positions=[[0.0, 0.0, 0.0]], cell=cell, pbc=True)
    Mn2 = Atoms("Mn2", scaled_positions=[[0.333, 0.333, 0.0], [0.6666, 0.666, 0.0]], cell=cell, pbc=True)

    return Au + Mn2

if __name__ == '__main__':

    alat = 3.1445
    factor = 2.5658

    npar = int(os.environ["SLURM_NNODES"])

    n_kpts_values = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 15, 19]

    primitive_unit_cell = make_primitive_Mn2Au_cell(alat, factor)

    if not os.path.exists("conv_files"):
        os.mkdir("conv_files")

    dataframe = pd.DataFrame()

    for n_kpts in n_kpts_values:

        print(f"Using {n_kpts} k points")

        VaspCalc = Vasp(npar=npar,
                        lorbit=11,
                        kpts=[n_kpts, n_kpts, n_kpts],
                        istart=0, # start from scratch
                        icharg=2, # default for istart=0
                        isif=2,
                        # nsim=2,
                        prec='Accurate',
                        encut=300,
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


        primitive_unit_cell.calc = VaspCalc
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

        dataframe = dataframe.append({"kpts" : n_kpts,
                                      "energy" : energy,
                                      "total_magmom" : primitive_unit_cell.get_magnetic_moment()
                                          }, ignore_index=True)

        dataframe.to_csv(f"kpts_conv_results.csv")
        primitive_unit_cell.write(f"conv_files/{n_kpts}_kpts_results.xyz")

        encut_values = [100, 200, 300, 400, 500, 600, 700, 800]

    dataframe = pd.DataFrame()
    for encut in encut_values:

        print(f"Using cutof {encut} eV")

        VaspCalc = Vasp(npar=npar,
                        lorbit=11,
                        kpts=[6, 6, 6],
                        istart=0, # start from scratch
                        icharg=2, # default for istart=0
                        isif=2,
                        nsim=2,
                        prec='Accurate',
                        encut=encut,
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


        primitive_unit_cell.calc = VaspCalc
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

        dataframe = dataframe.append({"encut" : encut,
                                      "energy" : energy,
                                      "total_magmom" : primitive_unit_cell.get_magnetic_moment()
                                          }, ignore_index=True)

        dataframe.to_csv(f"encut_conv_results.csv")
        primitive_unit_cell.write(f"conv_files/{encut}_encut_results.xyz")
