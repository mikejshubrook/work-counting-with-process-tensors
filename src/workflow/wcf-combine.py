
###################################################### combine.py ########################################################

# In this file we read through the WCF files that were created in wcf-run.py and stitch them together into a single csv file
# so that we can plot the WCF or fourier transform it to get the full work probability distribution.

#  The user runs this file by running the bash script "combine.sh".

# Standard library imports
import os

# Third-party imports
import numpy as np

# Path handling for reliable file saving/loading
CURRENT_FILE = os.path.abspath(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_FILE, "../../../"))  # Adjust based on your layout

DATA_DIR = os.path.join(PROJECT_ROOT, "data", "wcf-files")
WCF_DIR = os.path.join(PROJECT_ROOT, "data", "wcf-files-combined")

##### IMPORT PARAMETERS FROM BASH #####

# numerical parameters for iTEBD
PREC  = float(os.environ['PREC'])            # precision (threshold for SVD)
STEP_SIZE  = float(os.environ['STEP_SIZE'])  # step size for generalised time axis (trotter step)

# system hamiltonian parameters
EPSMAX = float(os.environ['EPSMAX'])         # maximum energy splitting
EPS0 = float(os.environ['EPS0'])             # minimum energy splitting
N=1
shift=True

# bath parameters
# extract alpha values from ALPHA_LIST
alpha_str = os.environ.get("ALPHA_LIST", "")
alpha_values = [float(x) for x in alpha_str.split()]
GAMMA = float(os.environ['GAMMA'])           # width of spectral density
W0 = float(os.environ['W0'])                 # peak location
BETA = float(os.environ['BETA'])             # inverse temperature
WC = float(os.environ['WC'])                 # cutoff frequency

# work counting parameters
S = int(os.environ['S']) # equilibration time - used for performing the first measurement in the two point measurement protocol
MSTART = int(os.environ['MSTART']) # starting value (in case you dont want to start at chi=0)
M = int(os.environ['M']) # number of steps along chi to take - what was the size of the job array when running wcf-run.sh   

# we convert the M values to chi, so that it is easy to compare WCFs with different values of STEP_SIZE
CHI0 = MSTART*STEP_SIZE
CHI_F =  (MSTART+M)*STEP_SIZE

# extract tau values from TAU_LIST
tp_str = os.environ.get("TP_LIST", "")
tp_values = [float(x) for x in tp_str.split()]

for ALPHA in alpha_values: # loop over coupling strengths
    for tau in tp_values:  # loop through protocol times
        for STA in [0, 1]:
            sta = (STA == 1)

            # Define folder name to grab files from
            folder = os.path.join(
                DATA_DIR,
                f"WCF_a{ALPHA}_G{GAMMA}_w{W0}_e{EPSMAX}_tp{tau}_sta{STA}_dt{STEP_SIZE}_p{PREC}_eq{S}",
                f"chi{MSTART * STEP_SIZE}"
            )
            # List of file names to extract data from
            file_names = [f"{folder}/wcf_m{m}.npy" for m in range(MSTART, MSTART+(M+1))] #TODO does this need the +1?

            # Output file where the combines data 
            output_file = f"{WCF_DIR}/WCF_a{ALPHA}_G{GAMMA}_w{W0}_e{EPSMAX}_t{tau}_sta{STA}_dt{STEP_SIZE}_p{PREC}_eq{S}_X{CHI0}_Xf{CHI_F}.txt"

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
                            # print(file_name, real, imag)
                            f_out.write(f"{real}\t{imag}\n")  # Write the complex number parts to the output file
                    else:
                        print(f"File {file_name} does not exist. Skipping.")

            print(f"All data has been successfully combined into {output_file}.")

