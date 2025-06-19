import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

if __name__ == "__main__":

    # Path handling for reliable file saving/loading
    CURRENT_FILE = os.path.abspath(__file__)
    PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_FILE, "../../"))  # Adjust based on your layout

    DATA_DIR = os.path.join(PROJECT_ROOT, "data", "dynamics-files")

    # numerical parameters for iTEBD
    PREC  = float(os.environ['PREC'])            # precision (threshold for SVD)
    STEP_SIZE  = float(os.environ['STEP_SIZE'])  # step size for generalised time axis (trotter step)

    # system hamiltonian parameters
    EPSMAX = float(os.environ['EPSMAX'])         # maximum energy splitting
    EPS0 = float(os.environ['EPS0'])             # minimum energy splitting

    # bath parameters
    GAMMA = float(os.environ['GAMMA'])           # width of spectral density
    W0 = float(os.environ['W0'])                 # peak location
    BETA = float(os.environ['BETA'])             # inverse temperature
    WC = float(os.environ['WC'])                 # cutoff frequency


    def plot_dynamics(a, tp, dt, p):
        """
        Plot the dynamics of the system based on the provided parameters.
        """
        
        fig, axs = plt.subplots(2, 3, figsize=(14, 6))
        
        for sta in [0, 1]:
            filename = f'{DATA_DIR}/dyns-a{a}_G10.0_w25.0_e25.0_tp{tp}_sta{sta}_dt{dt}_p{p}.csv'
            df = pd.read_csv(filename)
            
            # first row
            axs[0, 0].plot(df['tlist'], df['tr'], label=f'trace (sta={sta})')
            axs[0, 1].plot(df['tlist'], df['mag'], label=f'magnetization (sta={sta})')
            axs[0, 2].plot(df['tlist'], df['coh'], label=f'coherence (sta={sta})')
            # second row
            axs[1, 0].plot(df['tlist'], df['energy'], label=f'energy (sta={sta})')
            axs[1, 1].plot(df['tlist'], df['overlap_target'], label=f'overlap target (sta={sta})')
            axs[1, 2].plot(df['tlist'], df['overlap_instant'], label=f'overlap instant (sta={sta})')

        # formatting
        for ax in axs.flat:
            ax.set_xlabel('Time')
            ax.legend()
            ax.grid(alpha=0.2)
        
        plt.tight_layout()
        return

    alpha_str = os.environ.get("ALPHA_LIST", "")
    alpha_values = [float(x) for x in alpha_str.split()]

    tp_str = os.environ.get("TP_LIST", "")
    tp_values = [float(x) for x in tp_str.split()]

    alpha = alpha_values[0]  # Use the first alpha value for the plot
    tp = tp_values[0]        # Use the first tp value for the plot

    # Call the plotting function
    plot_dynamics(alpha, tp, STEP_SIZE, PREC)

    # Show the plot
    plt.show()