#!/bin/bash 

# Set number of threads (adjust to your CPU)
export OMP_NUM_THREADS=4

# counting parameter values
export MSTART=0 # starting value of counting integer
export M=50 # number of counting integers to calculate
export S=1 # time for equilibration

### Numerical parameters ###
export STEP_SIZE=0.01 # Trotter step size, must be tested for convergence
export PREC=7.0 # precision: 10^-{p} is the SVD threshold, must be tested for convergence

### Bath parameters ###
export BETA=1 # inverse temperature
export GAMMA=10.0 # spectral densit width
export W0=25.0 # spectral density peak location
export WC=5 # spectral density cutoff
export ALPHA=0.01 # coupling strength list, can be a single value or multiple values separated by spaces

### System parameters ###
export TP=1.0 # list of protocol times, can be a single value or multiple values separated by spaces
export EPSMAX=25.0 # maximum energy of the system
export EPS0=0.02 # minimum energy of the system

# Run scripts
python ../src/workflow/wcf.py
python ../src/workflow/wcf-combine.py
python wcf-test.py