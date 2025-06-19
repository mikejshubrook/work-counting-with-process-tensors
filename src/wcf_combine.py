# Standard library imports
import os 

# Third-party imports
import numpy as np  

# import numerical parameters from BASH

# numerical parameters for iTEBD
PREC  = float(os.environ['PREC'])            # precision (threshold for SVD)
STEP_SIZE  = float(os.environ['STEP_SIZE'])  # step size for generalised time axis (trotter step)

# system hamiltonian parameters
EPSTAU = float(os.environ['EPSTAU'])         # maximum energy splitting
EPS0 = float(os.environ['EPS0'])             # minimum energy splitting
N = int(os.environ['N'])                     # =1 for cyclic process, =2 for non-cyclic (Miller 2020)
SHIFT = int(os.environ['SHIFT'])             # =1 for shift in Hamiltonian, =0 for no shift
STA = int(os.environ['STA'])                 # =1 for shortcut to adiabaticity, =0 for no shortcut
TAU = float(os.environ['TAU'])               # protocol time

# bath parameters
ALPHA = float(os.environ['ALPHA'])           # system-bath coupling strength
GAMMA = float(os.environ['GAMMA'])           # width of spectral density
W0 = float(os.environ['W0'])                 # peak location
BETA = float(os.environ['BETA'])             # inverse temperature
WC = float(os.environ['WC'])                 # cutoff frequency

MSTART = int(os.environ['MSTART'])
M = int(os.environ['M'])
S = int(os.environ['S'])

if STA ==1:
    sta=True
else:
    sta=False


# Define folder name to grab files from
folder = f"wcf-files/WCF_a{ALPHA}_G{GAMMA}_w{W0}_e{EPSTAU}_tp{TAU}_SmoothSTA{STA}_dt{STEP_SIZE}_p{PREC}_eq{S}_N{N}_id{SHIFT}_g100_a0.1/chi{MSTART*STEP_SIZE}"
print(folder)
# List of file names to extract data from
file_names = [f"{folder}/wcf_m{m}.npy" for m in range(MSTART, MSTART+(M+1))]

# Output file where the combines data 
output_file = f"wcf-files/WCF_ALL_a{ALPHA}_G{GAMMA}_w{W0}_eps{EPSTAU}_tp{TAU}_SmoothSTA{STA}_dt{STEP_SIZE}_p{PREC}_eq{S}_N{N}_id{SHIFT}_g100_a0.1_X0{MSTART*STEP_SIZE}_Xf{(MSTART+M)*STEP_SIZE}.txt"
print(output_file)

print('here')

# Open the output file in write mode
with open(output_file, 'w') as f_out:
    print('file opened')
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


# combine with multiple starting values of chi
# used for convergence of WCF 


# for PREC in [9.0, 10.0]:

#     for chi_start in [0, 100, 500]:
#         MSTART = int(chi_start/STEP_SIZE)

#         # Define folder name to grab files from
#         folder = f"wcf-files/WCF_a{ALPHA}_G{GAMMA}_w{W0}_e{EPSTAU}_t{TAU}_dt{STEP_SIZE}_p{PREC}_eq{S}_N{N}_id{SHIFT}_STA{STA}/chi{chi_start}"

#         # List of file names to extract data from
#         file_names = [f"{folder}/wcf_m{m}.npy" for m in range(MSTART, MSTART+(M+1))]

#         # Output file where the combines data 
#         output_file = f"wcf-files/WCF_a{ALPHA}_G{GAMMA}_w{W0}_e{EPSTAU}_t{TAU}_dt{STEP_SIZE}_p{PREC}_eq{S}_N{N}_id{SHIFT}_STA{STA}_chi{chi_start}_{chi_start+M*STEP_SIZE}.txt"

#         # Open the output file in write mode
#         with open(output_file, 'w') as f_out:
#             # Loop through each file
#             for file_name in file_names:
#                 if os.path.exists(file_name):
#                     # Load the complex number from the .npy file
#                     complex_number = np.load(file_name)  # Assumes the file contains a single complex number

#                     # Check if it's indeed a complex number (optional but useful for error handling)
#                     if np.iscomplexobj(complex_number):
#                         real, imag = complex_number.real, complex_number.imag  # Extract real and imaginary parts
#                         f_out.write(f"{real}\t{imag}\n")  # Write the complex number parts to the output file

#         print(f"All data has been successfully combined into {output_file}.")
