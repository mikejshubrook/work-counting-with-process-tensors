import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

if __name__ == "__main__":

    # Path handling for reliable file saving/loading
    CURRENT_FILE = os.path.abspath(__file__)
    PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_FILE, "../../"))  # Adjust based on your layout

    DATA_DIR = os.path.join(PROJECT_ROOT, "data", "wcf-files-combined")

    # import numerical parameters from BASH

    # numerical parameters for iTEBD
    PREC  = float(os.environ['PREC'])            # precision (threshold for SVD)
    STEP_SIZE  = float(os.environ['STEP_SIZE'])  # step size for generalised time axis (trotter step)

    # system hamiltonian parameters
    EPSMAX = float(os.environ['EPSMAX'])         # maximum energy splitting
    EPS0 = float(os.environ['EPS0'])             # minimum energy splitting
    tp_str = os.environ.get("TP_LIST", "")
    tp_values = [float(x) for x in tp_str.split()]
    TP = tp_values[0]                            # time point for the WCF

    # bath parameters
    alpha_str = os.environ.get("ALPHA_LIST", "")
    alpha_values = [float(x) for x in alpha_str.split()]
    ALPHA= alpha_values[0]
    GAMMA = float(os.environ['GAMMA'])           # width of spectral density
    W0 = float(os.environ['W0'])                 # peak location
    BETA = float(os.environ['BETA'])             # inverse temperature
    WC = float(os.environ['WC'])                 # cutoff frequency

    # equilibration
    S = int(os.environ['S'])                     # equilibration time (same units as time)

    # counting
    MSTART = int(os.environ['MSTART']) # starting value (in case you dont want to start at chi=0)
    M = int(os.environ['M']) # number of steps along chi to take


    def plot_wcf(alpha, tp, dt, p, S, X0, Xf):
        """
        Plot the WCF (Work Correlation Function) based on the provided parameters.
        """
        fig, axs = plt.subplots(1, 2, figsize=(12, 4))
        
        for sta in [0, 1]:
            df = pd.read_csv(f'{DATA_DIR}/WCF_a{alpha}_G10.0_w25.0_e25.0_t{tp}_sta{sta}_dt{dt}_p{p}_eq{S}_X{X0}_Xf{Xf}.txt', sep='\t', header=None)
            chi_list = np.arange(0, len(df[0]))*dt

            # first row
            axs[0].plot(chi_list, df[0], label=f'Real (sta={sta})')

            # second row
            axs[1].plot(chi_list, df[1], label=f'Imag (sta={sta})')

        # formatting
        for ax in axs.flat:
            ax.set_xlabel('Counting parameter')
            ax.legend()
            ax.grid(alpha=0.2)
        
        plt.tight_layout()
        return

    plot_wcf(ALPHA, TP, STEP_SIZE, PREC, S, 0.0, float((M+MSTART)*STEP_SIZE))
    plt.show()