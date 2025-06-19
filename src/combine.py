
###################################################### combine.py ########################################################

# In this file we read through the WCF files that were created in wcf-run.py and stitch them together into a single csv file
# so that we can plot the WCF or fourier transform it to get the full work probability distribution.

#  The user runs this file by running the bash script "combine.sh".

# Standard library imports
import os

# Third-party imports
import numpy as np



##### IMPORT PARAMETERS FROM BASH #####

# numerical parameters for iTEBD
PREC  = float(os.environ['PREC'])            # precision (threshold for SVD)
STEP_SIZE  = float(os.environ['STEP_SIZE'])  # step size for generalised time axis (trotter step)

# system hamiltonian parameters
EPSMAX = float(os.environ['EPSMAX'])         # maximum energy splitting
EPS0 = float(os.environ['EPS0'])             # minimum energy splitting
TP = float(os.environ['TP'])                 # protocol time (tp)
STA = int(os.environ['STA'])                 # = 1 for STA, = 0 for no STA
N=1
shift=True

# bath parameters
ALPHA = float(os.environ['ALPHA'])           # coupling strength
GAMMA = float(os.environ['GAMMA'])           # width of spectral density
W0 = float(os.environ['W0'])                 # peak location
BETA = float(os.environ['BETA'])             # inverse temperature
WC = float(os.environ['WC'])                 # cutoff frequency

# work counting parameters
S = int(os.environ['S']) # equilibration time - used for performing the first measurement in the two point measurement protocol
MSTART = int(os.environ['MSTART']) # starting value (in case you dont want to start at chi=0)
M_TOTAL = int(os.environ['M_TOTAL']) # number of steps along chi to take - what was the size of the job array when running wcf-run.sh   

# we convert the M values to chi, so that it is easy to compare WCFs with different values of STEP_SIZE
CHI0 = MSTART*STEP_SIZE
CHI_F =  (MSTART+M_TOTAL)*STEP_SIZE

# Define folder name to grab files from
folder = f"wcf-files/WCF_a{ALPHA}_G{GAMMA}_w{W0}_e{EPSMAX}_tp{TP}_sta{STA}_dt{STEP_SIZE}_p{PREC}_eq{S}/chi{MSTART*STEP_SIZE}"

# List of file names to extract data from
file_names = [f"{folder}/wcf_m{m}.npy" for m in range(MSTART, MSTART+(M_TOTAL+1))]

# Output file where the combines data 
output_file = f"wcf-files/WCF_a{ALPHA}_G{GAMMA}_w{W0}_e{EPSMAX}_t{TP}_sta{STA}_dt{STEP_SIZE}_p{PREC}_eq{S}_X{CHI0}_Xf{CHI_F}.txt"

# Open the output file in write mode
with open(output_file, 'w') as f_out:
    # Loop through each file
    for file_name in file_names:
        if os.path.exists(file_name):
            
            # Load the complex number from the .npy file
            complex_number = np.load(file_name)  # Assumes the file contains a single complex number
            
            # Check if it's indeed a complex number (optional but useful for error handling)
            if np.iscomplexobj(complex_number):
                real, imag = complex_number.real, complex_number.imag  # Extract real and imaginary parts
                print(file_name, real, imag)
                f_out.write(f"{real}\t{imag}\n")  # Write the complex number parts to the output file
        else:
            print(f"File {file_name} does not exist. Skipping.")

print(f"All data has been successfully combined into {output_file}.")

