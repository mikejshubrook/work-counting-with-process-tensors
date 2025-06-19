
###################################################### moments.py ########################################################

# In this file we calculate work moments depending on user inputs from the bash script "moments.sh".
# The user can pass in a list of coupling strengths "ALPHA_LIST "and a list of protocol time values "TAU_LIST" 
# and this script will propagate the *generalised* system state, calculate the work characteristic function (WCF) at various points along
# the counting axis (chi), then use a finite difference order approximation to calculate derivatives of the WCF at chi=0, allowing for moments
# of the work probability distribution to be calculated.
# 
# The following moments are calculated:
# 1) mean
# 2) variance
# 3) skewness
# 4) kurtosis

# The user can also choose the order of the finite difference approximation (MAX_DIFF_ORDER). We make use of the findiff package (https://github.com/maroba/findiff)
# to calculate the finite difference coefficients for the derivatives. We store the moments for different orders of finite difference approximation so that
# we can check for convergence with the finite difference order.

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

DATA_DIR = os.path.join(PROJECT_ROOT, "data", "moments-files")
IF_DIR = os.path.join(PROJECT_ROOT, "data", "InfluenceFunctionals")



##### IMPORT PARAMETERS FROM BASH #####

# numerical parameters for iTEBD
PREC  = float(os.environ['PREC'])            # precision (threshold for SVD)
STEP_SIZE  = float(os.environ['STEP_SIZE'])  # step size for generalised time axis (trotter step)

# system hamiltonian parameters
EPSMAX = float(os.environ['EPSMAX'])         # maximum energy splitting
EPS0 = float(os.environ['EPS0'])             # minimum energy splitting
N=1 # cyclic process - used in sin function for the Hamiltonian
shift=True # shift the Hamiltonian to force ground state to be zero energy at all times

# bath parameters
GAMMA = float(os.environ['GAMMA'])           # width of spectral density
W0 = float(os.environ['W0'])                 # peak location
BETA = float(os.environ['BETA'])             # inverse temperature
WC = float(os.environ['WC'])                 # cutoff frequency

# work counting parameters
S = int(os.environ['S']) # equilibration time - used for performing the first measurement in the two point measurement protocol

# finite different approximation parameters
MAX_DIFF_ORDER = int(os.environ['MAX_DIFF_ORDER']) # maximum order of finite difference approximation
NTH_ORDER_DERIV = int(os.environ['NTH_ORDER_DERIV']) # order of the largest derivative we want to calculate (e.g. 1 for mean, 2 for mean and variance, etc.)

# extract alpha values from ALPHA_LIST
alpha_str = os.environ.get("ALPHA_LIST", "")
alpha_values = [float(x) for x in alpha_str.split()]
print("Alpha values:", alpha_values)

# extract tau values from TAU_LIST
tp_str = os.environ.get("TP_LIST", "")
tp_values = [float(x) for x in tp_str.split()]
print("tp values:", tp_values)

for ALPHA in alpha_values: # loop through coupling strengths

    # load the IF
    filename = os.path.join(IF_DIR, f'IF_a{ALPHA}_G{GAMMA}_w{W0}_dt{STEP_SIZE}_p{PREC}.pkl')
    with open(filename, 'rb') as file:
        itebd = pickle.load(file)
    
    for TP in tp_values: # loop through protocol times
        
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
            s = int(S/STEP_SIZE) # equilibration steps
            f = int(tp/STEP_SIZE) # physical time steps
            m_final = MAX_DIFF_ORDER + NTH_ORDER_DERIV # counting steps (allows us to calculate NTH_ORDER_DERIV up to MAX_DIFF_ORDER accuracy)

            # matrices of coefficients for finite difference approximation (assuming 4th order derivative is the largest we want)
            A1, A2, A3, A4 = get_deriv_coeffs(MAX_DIFF_ORDER, NTH_ORDER_DERIV)

            # empty array to store WCF values in
            wcf_chi = np.empty(m_final, dtype=complex)
            for m in range(m_final):  # loop through different steps along the counting axis
                # calculate WCF at final value of time (i.e. t=tp)
                wcf = process_iteration(itebd, s, f, m, Rho_0, STEP_SIZE, tp, eps0, eps_max, N, shift, sta)
                wcf_chi[m] = wcf

            # vectorized calculation of derivatives using finite difference coefficients from findiff
            d1 = (A1 @ wcf_chi)/STEP_SIZE
            d2 = (A2 @ wcf_chi)/STEP_SIZE**2
            d3 = (A3 @ wcf_chi)/STEP_SIZE**3
            d4 = (A4 @ wcf_chi)/STEP_SIZE**4

            # define moments of the work distrbution in terms of derivatives
            w1 = (-1j) * d1 # <W>
            w2 = (-1j)**2 * d2 # <W^{2}>
            w3 = (-1j)**3 * d3 # <W^{3}>
            w4 = (-1j)**4 * d4 # <W^{4}>

            # define mean and variance
            mean = w1
            variance = w2 - mean**2

            # avoid division by zero for skewness and kurtosis
            variance_safe = np.where(np.abs(variance) < 1e-12, np.nan, variance)

            # compute skewness and kurtosis with masking
            skewness_raw = (w3 - 3 * mean * variance - mean**3) / variance_safe**(3/2)
            kurtosis_raw = (w4 - 4 * mean * w3 - 6 * mean**2 * variance - 3 * mean**4) / variance_safe**2

            # replace NaNs with 0
            skewness = np.nan_to_num(skewness_raw, nan=0.0, posinf=0.0, neginf=0.0)
            kurtosis = np.nan_to_num(kurtosis_raw, nan=0.0, posinf=0.0, neginf=0.0)

            # Store all orders of finite difference approximation.
            # Store the real and imaginary parts too.
            # This will help us check for convergence with the finite difference order
            # The real part should converge to some valuel, and the imaginary parts should converge to zero.
            # store as dictionary
            data = {'mean_real':mean.real,
                    'variance_real': variance.real, 
                    'skewness_real': skewness.real,
                    'kurtosis_real': kurtosis.real,
                    'mean_imag': mean.imag,
                    'variance_imag': variance.imag,
                    'skewness_imag': skewness.imag,
                    'kurtosis_imag': kurtosis.imag,
                    }

            # Convert the dictionary to a pandas DataFrame
            df = pd.DataFrame(data)

            # Save the DataFrame to a CSV file
            output_path = os.path.join(DATA_DIR, f"moms-a{ALPHA}_G{GAMMA}_w{W0}_e{EPSMAX}_tp{TP}_sta{STA}_dt{STEP_SIZE}_p{PREC}_eq{S}.csv")

            # Ensure the output directory exists
            os.makedirs(os.path.dirname(output_path), exist_ok=True)

            # Save CSV
            df.to_csv(output_path, index=False)    

            print(f"Saved moments for ALPHA={ALPHA}, TP={TP}, STA={STA} to CSV.")