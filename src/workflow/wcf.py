# Standard library imports
import os 
import sys
import pickle  
import dill as pickle 

# Third-party imports
import numpy as np 
import pandas as pd 

# Get the directory of the current script (workflow)
current_dir = os.path.dirname(os.path.abspath(__file__))
# Get the parent directory (src)
src_dir = os.path.dirname(current_dir)
# Add src to sys.path
sys.path.append(src_dir)
# Local imports
from utils.LandauerErasure_functions import *

# Path handling for reliable file saving/loading
CURRENT_FILE = os.path.abspath(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_FILE, "../../../"))  # Adjust based on your layout

DATA_DIR = os.path.join(PROJECT_ROOT, "data", "wcf-files")
IF_DIR = os.path.join(PROJECT_ROOT, "data", "InfluenceFunctionals")
WCF_DIR = os.path.join(PROJECT_ROOT, "data", "wcf-files-combined")


# import numerical parameters from BASH

# numerical parameters for iTEBD
PREC  = float(os.environ['PREC'])            # precision (threshold for SVD)
STEP_SIZE  = float(os.environ['STEP_SIZE'])  # step size for generalised time axis (trotter step)

# system hamiltonian parameters
EPSMAX = float(os.environ['EPSMAX'])         # maximum energy splitting
EPS0 = float(os.environ['EPS0'])             # minimum energy splitting

# bath parameters
# extract alpha values from ALPHA_LIST
alpha_str = os.environ.get("ALPHA_LIST", "")
alpha_values = [float(x) for x in alpha_str.split()]
GAMMA = float(os.environ['GAMMA'])           # width of spectral density
W0 = float(os.environ['W0'])                 # peak location
BETA = float(os.environ['BETA'])             # inverse temperature
WC = float(os.environ['WC'])                 # cutoff frequency

# equilibration
S = int(os.environ['S'])                     # equilibration time (same units as time)

# counting
MSTART = int(os.environ['MSTART']) # starting value (in case you dont want to start at chi=0)
M = int(os.environ['M']) # number of steps along chi to take


# we convert the M values to chi, so that it is easy to compare WCFs with different values of STEP_SIZE
CHI0 = MSTART*STEP_SIZE
CHI_F =  (MSTART+M)*STEP_SIZE

N=1
shift=True # shift the Hamiltonian to force ground state to be zero energy at all times

# start at chi = 0 (remove this if you want to start elsewhere along chi axis)
m = MSTART+(M-1) # set counting integer

# scale parameters in terms of BETA
epstau = EPSMAX * BETA
eps0 = EPS0 * epstau
wc = WC * BETA

# extract tau values from TAU_LIST
tp_str = os.environ.get("TP_LIST", "")
tp_values = [float(x) for x in tp_str.split()]

# TODO: this code has been adapted from use on the HPC, which is why each value of m stores a separate
# WCF, which are then all combined at the end (due to using job arrays).
# It may be of use to instead remove the need for the combine script, and just save the WCFs in a single file.

for ALPHA in alpha_values: # loop over coupling strengths
    for tau in tp_values:  # loop through protocol times
        # for STA in [0, 1]:
        for STA in [0]:
            sta = (STA == 1)

            # define initial state (maximally mixed state)
            H = H_ERASURE(0, tau, eps0, epstau, N, shift, sta)
            evecs = np.linalg.eigh(H)[1]
            Rho_0 = 0.5 * (np.outer(evecs[:, 0], evecs[:, 0].conj()) + np.outer(evecs[:, 1], evecs[:, 1].conj()))

            # define steps along generalised time axis
            s = int(S/STEP_SIZE) # equilibration
            f = int(tau/STEP_SIZE) # time propagation

            # define filename of influence functional
            filename = os.path.join(IF_DIR, f'IF_a{ALPHA}_G{GAMMA}_w{W0}_dt{STEP_SIZE}_p{PREC}.pkl')

            # load IF
            with open(filename, 'rb') as file:
                itebd = pickle.load(file)

            print(f"Calculating WCF for ALPHA={ALPHA}, tp={tau}, STA={STA}, X0={MSTART*STEP_SIZE}, Xf={(MSTART+M)*STEP_SIZE}, S={S}")
            
            # --------------------------------------------------------------------------------------------------- #
            # if you are parallelising the code to run on multuple CPUs, uncomment the following section
            # you will then need to run 'wcf-combine.py' to combine the results from each CPU

            # for m in range(MSTART, MSTART + M):
            #     print(f"Calculating WCF for m={m} (chi={m*STEP_SIZE})")
            #     # calculate the WCF for this m
            #     wcf = process_iteration(itebd, s, f, m, Rho_0, STEP_SIZE, tau, eps0, epstau, N, shift, sta)
            #     # Construct folder path using os.path.join
            #     folder = os.path.join(
            #         DATA_DIR,
            #         f"WCF_a{ALPHA}_G{GAMMA}_w{W0}_e{EPSMAX}_tp{tau}_sta{STA}_dt{STEP_SIZE}_p{PREC}_eq{S}_p5",
            #         f"chi{MSTART*STEP_SIZE}")
            #     print(f"Saving WCF to {folder}")
            #     # Make sure the folder exists
            #     os.makedirs(folder, exist_ok=True)

            #     # Construct the full file path
            #     file_name = os.path.join(folder, f"wcf_m{m}.npy")

            #     # Save as .npy
            #     np.save(file_name, wcf)
            # --------------------------------------------------------------------------------------------------- #


            # --------------------------------------------------------------------------------------------------- #
            # if you are not parallelising the code, and want to save all WCFs in a single file, uncomment the following section
            # this will save the WCFs in a single file
            wcf_list = []  # list to store WCFs for each m
            for m in range(MSTART, MSTART + M):
                print(f"Calculating WCF for m={m} (chi={m*STEP_SIZE})")
                # calculate the WCF for this m
                wcf = process_iteration(itebd, s, f, m, Rho_0, STEP_SIZE, tau, eps0, epstau, N, shift, sta)
                # append the WCF to the list
                wcf_list.append(wcf)

            # Save the list of WCFs for this ALPHA and tau
            output_file = f"{WCF_DIR}/WCF_a{ALPHA}_G{GAMMA}_w{W0}_e{EPSMAX}_t{tau}_sta{STA}_dt{STEP_SIZE}_p{PREC}_eq{S}_p5_X{CHI0}_Xf{CHI_F}.txt"

            # Save the list of WCFs as text with real and imaginary parts in two columns
            with open(output_file, 'w') as f_out:
                for wcf in wcf_list:
                    f_out.write(f"{wcf.real:.16e}\t{wcf.imag:.16e}\n")
            # --------------------------------------------------------------------------------------------------- #


