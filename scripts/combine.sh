#!/bin/bash --login

#SBATCH -p multicore               # Partition name
#SBATCH -n 2                      # Number of cores
#SBATCH -t 7-0                     # Walltime (7 days)
#SBATCH -e ./output/%A_%a.err      # Standard error output
#SBATCH -o ./output/%A_%a.out      # Standard output

# number of cores
export OMP_NUM_THREADS=$SLURM_NTASKS

# counting parameter values
export M_TOTAL=20000 # number of counting parameters used (number of files to combine)
export MSTART=20000 # first value of M used in wcf job array
export S=5 # time for equilibration

### numerical parameters for the simulation ###
export STEP_SIZE=0.005 # step size for Trotter
export PREC=11.0 # SVD precision

### Physical parameters for the simulation ###
# system parameters
export EPSMAX=25.0 # maximum energy splitting
export EPS0=0.02 # minimum energy splitting
export TP=20.0  # protocol time
export STA=0 # =1 for STA, =0 for no STA (RAW)

# bath parameters
export ALPHA=0.16 # system-bath coupling strength
export GAMMA=10.0 # peak width
export W0=25.0 # peak position
export BETA=1 # inverse temperature
export WC=5 # cut-off frequency for spectral density (units of w0)

# load required modules
module load apps/binapps/anaconda3/2023.09

# RUN THIS FOR CONVERGENCE COMBINATION
python combine.py
