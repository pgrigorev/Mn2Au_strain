import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def linear(x, A):
    return A * x

def linear_fit(x, y, start=0, stop=None):

    popt, pcov = curve_fit(linear, x[start:stop], y[start:stop])

    print(popt)

    perr = np.sqrt(np.diag(pcov))
    print(perr)
    if stop is None:
        stop = -1
    x_new = np.arange(x[start], x[stop], 0.001)
    y_new = linear(x_new, *popt)

    return popt[0], x_new, y_new

if __name__ == '__main__':

    strain_results = pd.read_csv("../smaller_strain_110_results.csv", index_col=0)

    fig, ax = plt.subplots()

    ax.plot(strain_results.strain, strain_results.delta * 1.0e6, "o", label="$\Delta E = E(-110) - E(110)$")

    B, x_new, y_new = linear_fit(strain_results["strain"].values, strain_results.delta, start=1)
    ax.plot(x_new, y_new * 1.0e6, "C0--", label="fit with $\mathrm{\Delta E}$ = " + f"{B*1.0e3:.2f} x " + "$\mathrm{10^{-3}} \epsilon(110)$")

    ax.grid(True, linestyle="dashed", alpha=0.5)

    ax.legend()

    ax.set_xlabel("Strain along [110]")

    ax.set_ylabel("Energy ($\mathrm{\mu eV}$)")

    fig.savefig("smaller_strain_magnetostriction.pdf")
