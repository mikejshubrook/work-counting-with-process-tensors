###################################################### IF.py ########################################################

# In this file we create influence functional(s) depending on user inputs from a bash script.
# The user can pass in a list of alpha values (coupling strengths) and this python file will create the influence functional for each alpha value.
# The user also passes in numerical parameters (STEP_SIZE, PREC) and physical parameters (GAMMA, W0, BETA, WC) via the bash script.

# Standard library imports
import os
import sys
import pickle
import dill as pickle

# Third-party imports
import numpy as np

# Local imports
current_dir = os.path.dirname(os.path.abspath(__file__)) # Get the directory of the current script (workflow)
src_dir = os.path.dirname(current_dir) # Get the parent directory (src)
sys.path.append(src_dir) # Add src to sys.path
from utils.LandauerErasure_functions import *

# Get absolute path to project root (2 levels up from this file)
CURRENT_FILE = os.path.abspath(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_FILE, "../../../"))

# Convenience: path to data folder
DATA_DIR = os.path.join(PROJECT_ROOT, "data", "InfluenceFunctionals")

# import numerical parameters from BASH
STEP_SIZE  = float(os.environ['STEP_SIZE'])  # step size for generalised time axis (trotter step)
PREC  = float(os.environ['PREC'])            # precision (threshold for SVD)kjz     

# import physical parameters from BASH
GAMMA = float(os.environ['GAMMA'])           # width of spectral density
W0 = float(os.environ['W0'])                 # peak location
BETA = float(os.environ['BETA'])             # inverse temperature
WC = float(os.environ['WC'])                 # cutoff frequency

# extract alpha values from ALPHA_LIST
alpha_str = os.environ.get("ALPHA_LIST", "")
alpha_values = [float(x) for x in alpha_str.split()]
print("Alpha values:", alpha_values)

# define parameters for the influence functional (which are independent of alpha)
# TODO: it might be better to define these in the bash script so they are easily adjustablekk
hard_freq_cutoff = W0*WC                # hard cutoff for frequency integrals
T_max = 5                               # bath memory time
d = 3                                   # regularizing function gradient
tr = 2.5                                # regularizing function location
MEMORY_STEPS = int(T_max / STEP_SIZE)   # memory steps

# loop over alpha values and create influence functional in each case
for ALPHA in alpha_values:
                
    # define filename of influence functional that we will create
    filename = f'data/InfluenceFunctionals/IF_a{ALPHA}_G{GAMMA}_w{W0}_dt{STEP_SIZE}_p{PREC}.pkl'

    # define the bath spectral density for this value of alpha
    J = lambda w:  (ALPHA*GAMMA*W0**2 *w) /( (W0**2 - w**2)**2 + (GAMMA*w)**2 )


    filename = os.path.join(DATA_DIR, f"IF_a{ALPHA}_G{GAMMA}_w{W0}_dt{STEP_SIZE}_p{PREC}.pkl")

    # Ensure directory exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    # # Check if the file already exists to avoid remaking and overwriting
    if os.path.exists(filename):
        print(f'{filename} already exists. Skipping creation of influence functional.')

    else:
        # create the influence functional
        itebd = find_IF(PREC, STEP_SIZE, 
                        MEMORY_STEPS, T_max, 
                        hard_freq_cutoff, J, BETA, 
                        d, tr, check_regularisation=False)

        with open(filename, 'wb') as f:
            pickle.dump(itebd, f)

        print('Created the influence function and saved to: ', filename)