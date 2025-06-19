#!/bin/bash

# Set number of threads (adjust to your CPU)
export OMP_NUM_THREADS=4

### Numerical parameters ###
export STEP_SIZE=0.01
export PREC=7.0

### Physical parameters ###
export GAMMA=10.0
export W0=25.0
export BETA=1
export WC=5
export ALPHA_LIST="0.02"
export TP_LIST="1 2 3"
export EPSMAX=25.0
export EPS0=0.02

### Parameters for moments calculations ###
export S=1
export NTH_ORDER_DERIV=4
export MAX_DIFF_ORDER=4

# Run scripts
python ../src/IF.py
python ../src/dynamics.py
python ../src/moments.py
