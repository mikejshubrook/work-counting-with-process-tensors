#!/bin/bash --login

#SBATCH -p multicore               # Partition name
#SBATCH -n 2                      # Number of cores
#SBATCH -t 7-0                     # Walltime (7 days)
#SBATCH -a 1-2                     # Job array range
#SBATCH -e ./output/%A_%a.err      # Standard error output
#SBATCH -o ./output/%A_%a.out      # Standard output
#SBATCH --mail-type=All            # Mail events
#SBATCH --mail-user=mike.shubrook@manchester.ac.uk

# Number of cores
export OMP_NUM_THREADS=$SLURM_NTASKS

# if you dont want to run the dynamics.py file, comment out the line "python dynamics.py"
# if you dont want to run the moments.py file, comment out the line "python moments.py"

### Numerical parameters for the simulations ###
export STEP_SIZE=0.01 # step size for Trotter
export PREC=$(echo "scale=2; 10.0 + (1 * ($SLURM_ARRAY_TASK_ID-1))" | bc) # SVD precision 

### Physical parameters for the simulations ###
# bath parameters
export GAMMA=10.0 # peak width
export W0=25.0 # peak position
export BETA=1 # inverse temperature
export WC=5 # cut-off frequency for spectral density (units of w0)
# export ALPHA_LIST="0.01 0.04 0.16"
# export ALPHA_LIST="0.02 0.08 0.32"
export ALPHA_LIST="0.32"

# system parameters (ignored by IF.py)
export EPSMAX=25.0 # maximum energy splitting
export EPS0=0.02 # minimum energy splitting
# export TP_LIST="0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95 1.05 1.15 1.25 1.35 1.45 1.55 1.65 1.75 1.85 1.95"
# export TP_LIST="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.2 5.4 5.6 5.8 6.0 6.2 6.4 6.6 6.8 7.0 7.2 7.4 7.6 7.8 8.0 8.2 8.4 8.6 8.8 9.0 9.2 9.4 9.6 9.8 10.0 10.2 10.4 10.6 10.8 11.0 11.2 11.4 11.6 11.8 12.0 12.2 12.4 12.6 12.8 13.0 13.2 13.4 13.6 13.8 14.0 14.2 14.4 14.6 14.8 15.0 15.2 15.4 15.6 15.8 16.0 16.2 16.4 16.6 16.8 17.0 17.2 17.4 17.6 17.8 18.0 18.2 18.4 18.6 18.8 19.0 19.2 19.4 19.6 19.8 20.0"
# export TP_LIST="0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.2 5.4 5.6 5.8 6.0 6.2 6.4 6.6 6.8 7.0 7.2 7.4 7.6 7.8 8.0 8.2 8.4 8.6 8.8 9.0 9.2 9.4 9.6 9.8 10.0 10.2 10.4 10.6 10.8 11.0 11.2 11.4 11.6 11.8 12.0 12.2 12.4 12.6 12.8 13.0 13.2 13.4 13.6 13.8 14.0 14.2 14.4 14.6 14.8 15.0 15.2 15.4 15.6 15.8 16.0 16.2 16.4 16.6 16.8 17.0 17.2 17.4 17.6 17.8 18.0 18.2 18.4 18.6 18.8 19.0 19.2 19.4 19.6 19.8 20.0"
export TP_LIST="1.0"
# parameters needed for moments calculations (ignored by IF.py and dynamics.py)
export S=5 # time for equilibration
export NTH_ORDER_DERIV=4 # nth order derivative for finite difference approximation (up to and including 4th order i.e. kurtosis
export MAX_DIFF_ORDER=10 # mth order accuracy of finite difference approximation

# load required modules
module load apps/binapps/anaconda3/2023.09

# run the influence functional code
# python IF.py

# run the dynamics file
python dynamics-states.py

# run the moments file
# python moments.py