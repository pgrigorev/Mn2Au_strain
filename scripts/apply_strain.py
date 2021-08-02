import os

import pandas as pd

from ase import Atoms
from ase.calculators.vasp import Vasp

if __name__ == '__main__':

    primitive_unit_cell = read("strain_filtered_unit_cell.xyz")
    npar = int(os.environ["SLURM_NNODES"])


    moments = np.array([[0., 0., 0.],
                        [0., 0., -3.5],
                        [0., 0., 3.5]])

    primitive_unit_cell.set_initial_magnetic_moments(moments)


    strains = [-0.05, -0.01, -5e-3, -1e-3, 0.0, 1e-3, 5e-3, 0.01, 0.05]

    results = pd.DataFrame()

    configs_folder = "configurations"

    if not os.path.exists(configs_folder):
        os.mkdir(configs_folder)

    for strain in strains:

        configuration = primitive_unit_cell.copy()
        configutation.cell[2][0] += strain * primitive_unit_cell[2][0]
        configutation.cell[2][1] += strain * primitive_unit_cell[2][1]

        VaspCalc = Vasp(npar=npar,
                        lorbit=11,
                        kpts=[11, 11, 11],
                        istart=1, # read collinear wavecar
                        # icharg=0, # non self consistent
                        isif=2,
                        # nsim=2,
                        lmaxmix=4,
                        lwave=False, # do not write wavecar
                        lcharg=False, # do not write CHGCAR
                        nbands=2 * 16,
                        prec='Accurate',
                        encut=600,
                        ediff=1.e-7,
                        lreal=False,
                        nelm=200,
                        nelmin=5,
                        algo="Normal",
                        isym=-1,
                        ismear=1,
                        sigma=0.1,
                        ispin=2, # spin polarised calculation
                        lsorbit=True, # spin orbit coupling
                        saxis=[1, 1, 0],
                        xc='PBE')

        configutation.calc = VaspCalc

        energy = configutation.get_potential_energy()


        results = resluts.append({"strain": strain,
                                  "energy": energy}, ignore_index=True)
        results.to_csv("strain_results.csv")

        primitive_unit_cell.wrte(configs_folder +
                                f"/config_strain{strain}.xyz")
