
###################################################### dynamics.py ########################################################

# In this file we run the system dynamics depending on user inputs from the bash script, after creating the influence functional in IF.py.
# The user can pass in a list of coupling strengths "ALPHA_LIST "and a list of protocol time values "TP_LIST" 
# and this script will propagate the system state and save the following properties to file:
# 1) trace of the density operator
# 2) magnetization
# 3) coherence
# 4) expectation of system Hamiltonian
# 5) Instantaneous overlap with the time dependent Hamiltonian
# 6) Overlap with the ground state of the final Hamiltonian - useful for erasure protocols to check for erasure 'success'.

# The user also passes in numerical parameters (STEP_SIZE, PREC) and physical parameters (GAMMA, W0, BETA, WC) via the bash script.

# It is recommended to use a job array to run this script for diferent PREC values to ensure convergence of the results wrt PREC,
# and then to repeat this for different STEP_SIZE values to ensure convergence of the results wrt STEP_SIZE.

# Standard library imports
import os
import pickle
import dill as pickle

# Third-party imports
import numpy as np
import pandas as pd

# Local imports
from utils.LandauerErasure_functions import *


# Path handling for reliable file saving/loading
CURRENT_FILE = os.path.abspath(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_FILE, "../../"))  # Adjust based on your layout

DATA_DIR = os.path.join(PROJECT_ROOT, "data", "dynamics-files")
IF_DIR = os.path.join(PROJECT_ROOT, "data", "InfluenceFunctionals")


##### IMPORT PARAMETERS FROM BASH #####

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

# extract alpha values from ALPHA_LIST
alpha_str = os.environ.get("ALPHA_LIST", "")
alpha_values = [float(x) for x in alpha_str.split()]
print("Alpha values:", alpha_values)

# extract tau values from TAU_LIST
tp_str = os.environ.get("TP_LIST", "")
tp_values = [float(x) for x in tp_str.split()]
print("tp values:", tp_values)

shift = True # shift the Hamiltonian to force ground state to be zero energy at all times
N=1 # cyclic process - used in sin function for the Hamiltonian

# loop over all alpha values
for ALPHA in alpha_values:
    
    # define filename of influence functional (IF)
    # filename = f'InfluenceFunctionals/IF_a{ALPHA}_G{GAMMA}_w{W0}_dt{STEP_SIZE}_p{PREC}.pkl'
    filename = os.path.join(IF_DIR, f'IF_a{ALPHA}_G{GAMMA}_w{W0}_dt{STEP_SIZE}_p{PREC}.pkl')

    # load the influence functional
    with open(filename, 'rb') as file:
        itebd = pickle.load(file)
    
    # for each alpha value, loop over all tau values
    for TP in tp_values:
        
        # scale physical parameters in terms of BETA
        eps_max = EPSMAX * BETA
        eps0 = EPS0 * eps_max
        wc = WC * BETA
        tp = TP

        # for each alpha and tau combination, we want to check the dynamics for both STA and non-STA
        for STA, sta in [(1, True), (0, False)]: # STA or raw Hamiltonian (with labels for saving file names)

            # define initial state of TLS (maximally mixed state of the initial Hamiltonian eigenstates)
            H = H_ERASURE(0, tp, eps0, eps_max, N, shift, sta)
            evecs = np.linalg.eigh(H)[1]
            Rho_0 = 0.5 * (np.outer(evecs[:, 0], evecs[:, 0].conj()) + np.outer(evecs[:, 1], evecs[:, 1].conj()))

            # define the number of steps along the generalised time axis (tau)
            f = int(tp/STEP_SIZE) # physical time steps

            # define the time list
            tlist = np.array(np.arange((f+1)) * STEP_SIZE)

            # calculate the system dynamics
            states = system_dynamics(itebd, 0, f, 0, Rho_0, STEP_SIZE, tp, eps0, eps_max, N, shift, sta)

            # calculate trace, magnetization, coherence
            tr, mag, coh = trace_mag_coh(states)

            # calculate system hamiltonian expectation value
            energy = system_energy(states, tlist, tp, eps0, eps_max, N, shift, sta)

            # calcualte overlaps (target and instantaneous)
            overlap_target, overlap_instant = overlap(states, tlist, tp, eps0, eps_max, N, shift, sta)

            # store results in a dictionary
            results = {
                'tlist': tlist.real,
                'tr': tr.real,
                'mag': mag.real,
                'coh': coh.real,
                'energy': energy.real,
                'overlap_target': overlap_target.real,
                'overlap_instant': overlap_instant.real
            }

            # Convert dict of resuilts to a pandas DataFrame
            df = pd.DataFrame(results)

            # Save DataFrame to a CSV file
            # df.to_csv(f"data/dynamics-files/dyns-a{ALPHA}_G{GAMMA}_w{W0}_e{EPSMAX}_tp{TP}_sta{STA}_dt{STEP_SIZE}_p{PREC}.csv", index=False)
            output_path = os.path.join(DATA_DIR, f"dyns-a{ALPHA}_G{GAMMA}_w{W0}_e{EPSMAX}_tp{TP}_sta{STA}_dt{STEP_SIZE}_p{PREC}.csv")

            # Ensure the output directory exists
            os.makedirs(os.path.dirname(output_path), exist_ok=True)

            # Save CSV
            df.to_csv(output_path, index=False)
            
            print(f"Saved dynamics for ALPHA={ALPHA}, TP={TP}, STA={STA} to CSV.")