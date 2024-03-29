import os

import pandas as pd

from ase import Atoms
from ase.calculators.vasp import Vasp

def get_angle(a, b):
    a = np.array(a)
    b = np.array(b)
    unit_a = a / np.linalg.norm(a)
    unit_b = b / np.linalg.norm(b)
    dot_product = np.dot(unit_a, unit_b)
    return np.degrees(np.arccos(dot_product))

if __name__ == '__main__':

    primitive_unit_cell.read("strain_filtered_unit_cell.xyz")
    npar = int(os.environ["SLURM_NNODES"])


    moments = np.array([[0., 0., 0.],
                        [0., 0., -3.5],
                        [0., 0., 3.5]])

    primitive_unit_cell.set_initial_magnetic_moments(moments)


    directions = [[-1, 1, 0],
                   [-9, 10, 0],
                   [3, 4, 0],
                   [-2, 3,0],
                   [-1, 2, 0],
                   [-1, 3,0],
                   [-1, 4, 0],
                   [-1, 10, 0],
                   [-1, 20, 0],
                   [0, 1 ,0],
                   [1, 20,0],
                   [1, 10, 0],
                   [1, 4, 0],
                   [1, 3, 0],
                   [1, 2, 0],
                   [2, 3, 0],
                   [3, 4, 0],
                   [9, 10, 0],
                   [1, 1, 0]]

    results = pd.DataFrame()

    configs_folder = "configurations"

    if not os.path.exists(configs_folder):
        os.mkdir(configs_folder)

    for direction in directions:

        VaspCalc = Vasp(npar=npar,
                        lorbit=11,
                        kpts=[11, 11, 11],
                        istart=1, # read collinear wavecar
                        icharg=11, # non self consistent
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
                        saxis=direction,
                        xc='PBE')

        primitive_unit_cell.calc = VaspCalc

        energy = primitive_unit_cell.get_potential_energy()
        angle = get_angle(direction, [1, 1, 0])
        print(f"Angle from [110]: {angle}")
        print(f"Energy: {energy}")

        results = resluts.append({"angle": angle,
                                  "tetta": get_angle(direction, [1,0,0]),
                                  "energy": energy,
                                  "dir_x": direction[0],
                                  "dir_y": direction[1],
                                  "dir_z": direction[2]}, ignore_index=True)
        results.to_csv("110_anisotropy_results.csv")

        primitive_unit_cell.wrte(configs_folder +
                                f"/config_direction{str(direction).replace(", ", "_")}.xyz")
