# Work counting with process tensors

This project calculates numerically exact quantum work statistics with process tensors.

# Project Overview
This code was used to calculate the results in (link to paper).


# Initial setup
To run the code, follow these instructions:
TODO: basic setup, requirements, etc

Create the data folders by running (within 'data' folder) 
"bash setup_data_folder.sh"

# Tests
Run the dynamics test by running (within 'tests' folder) 
"bash run-dynamics-tests.sh" 
It will first create an influence functional and then run the system dynamics, and finally it will plot these dynamics.

Run the WCF test by running (within 'tests' folder) 
"bash run-wcf-tests.sh" 
This will calculate the WCF using the same parameters used to calculate the influence functional and system dynamics from the previous test. 

To see if your tests have run properly, compare the outputted figures to the files within the 'tests' folder names 'expected-dynamics.png' and 'expected-wcf.png'

TODO: add in moments test as well.

# List of parameters
This project has two main categories of parameters that the user should know about, being physical parameters and numerical parameters. 

The physical parameters tell us about the system that we are trying to simulate, e.g. energy scales, time, etc.
Numerical parameters are what we use to simulate a given physical situation, e.g. time step, SVD threshold

The goal is to get 'convergence' of the results with repsect to the numerical paraemters, for a given set of physical parameters. When this has been acheived, 
the results are converged and therefore the results that we are seeing are not due to the way we are simulating this system, but are due to the results of the 
simulation itself.

Below are a breakdown of the physical and numerical parameters used in this project. I will cover the parameters which are used within the bash scripts, 
and leave the description of other parameters used elsewhere to the files in which they are used.

Physical parameters: broken down into system, bath, and shared parameters
System: 
EPSMAX: (float) maximum energy splitting of the TLS
EPS0: (float) minimum energy splitting of the TLS
TP: (float) protocol time, i.e. how long to drive the TLS for

Bath: 
(note, here we use an underdamped drude-lorentz spectral density, it is possible to change this by changing the spectral density function within src/LandauerErasure_functions.py and the corresponding functions which depend on this)
BETA: (float) inverse temperature, and the base unit (i.e. all other parameters are measured in units of BETA)
GAMMA: (float) width of the spectral density
W0: (float) peak location of the spectral density
WC: (float) cutoff for the spectral density
ALPHA: (float) coupling strength between the system and bath

Numerical parameters: (i.e. simulation parameters)
STEP_SIZE: Trotter step size, used for discretization of the generalised time axis
PREC: threshold of SVD, in units of 10^{-PREC}

Other parameters:
S: (float) equilibration time, used to perform the first measurement within the TPMP, must be chosen to be large enough to have converged results
M: (int) counting integer, used to calcualte the WCF along the chi axis, in units of M*STEP_SIZE
MAX_DIFF_ORDER: (int) used when calculating moments through the finite difference approximation, to calculate those approximations to higher order accuracy (should be an integer and multiple of 2)

# Example usage
I have uploaded video tutorials on my YouTube channel explaining both how the code works, and how to use it. 
Please use the following link to see the playlist of videos associated with this project:
(link to playlist)


