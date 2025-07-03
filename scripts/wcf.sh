#!/bin/bash 

# This script allows the user to calculate the work characteristic function (WCF) for a given set of physical and numerical parameters.
# First, the influence functional is calculated, then the work characteristic function data points, and finally these data points combined.
# This is an artefact of the worflow, that was developed to run on HPC in parallel.


# Set number of threads (adjust to your CPU)
export OMP_NUM_THREADS=11

# counting parameter values
export MSTART=0 # starting value of counting integer
export S=5 # time for equilibration

### Numerical parameters ###
export STEP_SIZE=0.01 # Trotter step size, must be tested for convergence
export PREC=7.0 # precision: 10^-{p} is the SVD threshold, must be tested for convergence

# Choose maximum value of the counting paramter chi
X=200.0

# Calculate M = X / STEP_SIZE, as an integer
export M=$(printf "%.0f" $(echo "$X / $STEP_SIZE" | bc -l)) # number of counting integers to calculate

### Bath parameters ###
export BETA=1 # inverse temperature
export GAMMA=10.0 # spectral densit width
export W0=25.0 # spectral density peak location
export WC=5 # spectral density cutoff
export ALPHA_LIST="0.02" # coupling strength list, can be a single value or multiple values separated by spaces

### System parameters ###
export TP_LIST="2.0" # list of protocol times, can be a single value or multiple values separated by spaces
export EPSMAX=25.0 # maximum energy of the system
export EPS0=0.02 # minimum energy of the system

# Run py files
python ../src/workflow/IF.py
python ../src/workflow/wcf.py
python ../src/workflow/wcf-combine.py



