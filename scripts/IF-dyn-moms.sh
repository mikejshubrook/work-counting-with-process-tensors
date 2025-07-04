#!/bin/bash

# This script allows the user to calculate system dynamics and work moments for a given set of physical and numerical parameters.
# First, the influence functional is calculated, then the system dynamics, and finally the work moments.

# Set number of threads (adjust to your CPU)
export OMP_NUM_THREADS=4

### Numerical parameters ###
export STEP_SIZE=0.01 # Trotter step size, must be tested for convergence
export PREC=7.0 # precision: 10^-{p} is the SVD threshold, must be tested for convergence

### Bath parameters - required for IF.py, dynamics.py and moments.py ###
export BETA=1 # inverse temperature
export GAMMA=10.0 # spectral densit width
export W0=25.0 # spectral density peak location
export WC=5 # spectral density cutoff
export ALPHA_LIST="0.02" # coupling strength list, can be a single value or multiple values separated by spaces

### System parameters - only required for dynamics.py and moments.py ###
export TP_LIST="2.0 20.0" # list of protocol times, can be a single value or multiple values separated by spaces
export EPSMAX=25.0 # maximum energy of the system
export EPS0=0.02 # minimum energy of the system

### Moments parameters, only required for moments.py ###
export S=1 # equilibration time, must be tested for convergence
export MAX_DIFF_ORDER=12 # maximum order of the finite difference approximation for derivatives

# Run py files
python ../src/workflow/IF.py
python ../src/workflow/dynamics.py
# python ../src/workflow/moments.py
