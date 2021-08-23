import os
import numpy as np

from shutil import copyfile
import pandas as pd

from ase.io import read
from ase.calculators.vasp import Vasp

if __name__ == '__main__':

    primitive_unit_cell = read("strain_filtered_unit_cell.xyz")
    npar = int(os.environ["SLURM_NNODES"])


    moments = np.array([[0., 0., 0.],
                        [0., 0., -3.5],
                        [0., 0., 3.5]])

    primitive_unit_cell.set_initial_magnetic_moments(moments)


    strains = [-0.02, -0.015,  -0.01, -5e-3,
               0.0, 
               5e-3, 0.01, 0.015, 0.02]
    if os.path.exists("strain_results.csv"):
        results = pd.read_csv("strain_results.csv", index_col=0)
    else:
        results = pd.DataFrame()

    configs_folder = "configurations"

    if not os.path.exists(configs_folder):
        os.mkdir(configs_folder)

    starting_index = len(results)
    print(f"Starging from strained configuration {starting_index}")
    for strain in strains[starting_index:]:

        configuration = primitive_unit_cell.copy()
        cell = configuration.cell.copy()
        cell[0][0] += strain * primitive_unit_cell.cell[0][0] / np.sqrt(2.0)
        cell[0][1] += strain * primitive_unit_cell.cell[0][1] / np.sqrt(2.0)
        cell[1][0] += strain * primitive_unit_cell.cell[1][0] / np.sqrt(2.0)
        cell[1][1] += strain * primitive_unit_cell.cell[1][1] / np.sqrt(2.0)
        configuration.set_cell(cell, scale_atoms=False)
        
        VaspSTD = Vasp(npar=npar,
                        kpar=10,
                        lorbit=11,
                        kpts=[25, 25, 25],
                        istart=0, # start non-collinear run
                        # icharg=11, # non self consistent
                        isif=2,
                        # nsim=2,
                        lmaxmix=4,
                        # lwave=False, # do not write wavecar
                        # lcharg=False, # do not write CHGCAR
                        # nbands=2 * 20,
                        prec='Accurate',
                        encut=600,
                        ediff=1.e-9,
                        lreal=False,
                        nelm=200,
                        nelmin=5,
                        algo="Normal",
                        isym=-1,
                        ismear=-5,
                        # sigma=0.1,
                        ispin=2, # spin polarised calculation
                        # lsorbit=True, # spin orbit coupling
                        addgrid=True,
                        gga_compat=False,
                        # saxis=[-1, 1, 0],
                        magmom=[0.0, -3.5, 3.5],
                        xc='PBE')

        configuration.calc = VaspSTD
        os.environ["VASP_COMMAND"] = "srun vasp_std"
        collinear_energy = configuration.get_potential_energy() # generate WAVECAR, CHG, CHGCAR
        print(f"collinear_energy: {collinear_energy}")
        os.environ["VASP_COMMAND"] = "srun vasp_ncl"

        VaspCalc = Vasp(npar=npar,
                        kpar=10,
                        lorbit=11,
                        kpts=[25, 25, 25],
                        istart=1, # read collinear wavecar
                        icharg=11, # non self consistent
                        isif=2,
                        # nsim=2,
                        lmaxmix=4,
                        lwave=False, # do not write wavecar
                        lcharg=False, # do not write CHGCAR
                        nbands=2 * 20,
                        prec='Accurate',
                        encut=600,
                        ediff=1.e-9,
                        lreal=False,
                        nelm=200,
                        nelmin=5,
                        algo="Normal",
                        isym=-1,
                        ismear=-5,
                        # sigma=0.1,
                        ispin=2, # spin polarised calculation
                        lsorbit=True, # spin orbit coupling
                        addgrid=True,
                        gga_compat=False,
                        saxis=[1, 1, 0],
                        xc='PBE')

        configuration.calc = VaspCalc

        perpendicular_energy = configuration.get_potential_energy()

        copyfile("OUTCAR", configs_folder +  f"/OUTCAR_strain_{strain}_perpencular")
        VaspCalc = Vasp(npar=npar,
                        kpar=10,
                        lorbit=11,
                        kpts=[25, 25, 25],
                        istart=1, # read collinear wavecar
                        icharg=11, # non self consistent
                        isif=2,
                        # nsim=2,
                        lmaxmix=4,
                        lwave=False, # do not write wavecar
                        lcharg=False, # do not write CHGCAR
                        nbands=2 * 20,
                        prec='Accurate',
                        encut=600,
                        ediff=1.e-9,
                        lreal=False,
                        nelm=200,
                        nelmin=5,
                        algo="Normal",
                        isym=-1,
                        ismear=-5,
                        # sigma=0.1,
                        ispin=2, # spin polarised calculation
                        lsorbit=True, # spin orbit coupling
                        addgrid=True,
                        gga_compat=False,
                        saxis=[-1, 1, 0],
                        xc='PBE')

        configuration.calc = VaspCalc

        parallel_energy = configuration.get_potential_energy()
        delta = perpendicular_energy - parallel_energy
        results = results.append({"strain": strain,
                                  "collinear_energy": collinear_energy,
                                  "perpendicular_energy": perpendicular_energy,
                                  "parallel_energy": parallel_energy,
                                  "delta": delta}, ignore_index=True)

        results.to_csv("strain_results.csv")

        # configuration.write(configs_folder +
        #                        f"/config_strain_{strain}_.xyz")
        copyfile("vasp.out", configs_folder +  f"/vasp_strain_{strain}_.out")
        copyfile("OUTCAR", configs_folder +  f"/OUTCAR_strain_{strain}_parallel")
        copyfile("POSCAR", configs_folder + f"/POSCAR_strain_{strain}")
