
# Third-party imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
import findiff

# Local imports
from utils.WC_TEMPO import time_dependent_TEMPO
from utils.bath_functions import BathCorrelation
from utils.operators import sz, sy, sx, spre, spost


#################### CREATE INFLUENCE FUNCTIONAL WITH ITEBD ####################

def find_IF(prec, STEP_SIZE, MEMORY_STEPS, T_max, wc,  J, BETA, delta, tr, check_regularisation = False):

    """
    Create the influence functional using the ITEBD method.

    Parameters:
        prec (float): Precision for the computation.
        STEP_SIZE (float): Step size for the time evolution, used as the Trotter step size.
        MEMORY_STEPS (int): Number of memory steps.
        T_max (float): Maximum time for the bath correlation function.
        wc (float): Cutoff frequency for the bath correlation function.
        J (function): Spectral density function.
        BETA (float): Inverse temperature.
        delta (float): Regularization parameter.
        tr (float): Regularization time.
        check_regularisation (bool): If True, check the regularization of the bath correlation function.
        
    Returns:
        itebd (object): The influence functional object.
    """

    # Define the precision for the computation
    precision = 10**(-prec)

    # define the IF object
    itebd = time_dependent_TEMPO(np.diagonal(sz()), STEP_SIZE, MEMORY_STEPS)
    BCF = BathCorrelation()

    #Define the regularizing function:
    taper_func = lambda t: 1/(1 + np.exp(delta * (t - tr)))

    # Define the regularization parameters to be passed to the bath correlation function
    regularize = [True, taper_func]

    # Construct the bath correlation function
    BCF.construct_bath_correlation(J, BETA, T_max, MEMORY_STEPS, 20 * wc, regularize)

    # if the user wants to check the regularization, plot the bath correlation functions with 
    # and without regularization 
    if check_regularisation:
        
        # create a bath correlation object without regularization
        BCF_NOREG = BathCorrelation()
        BCF_NOREG.construct_bath_correlation(J, BETA, T_max, MEMORY_STEPS, 20 * wc, regularize)
        
        # plot the bath correlation function WITH and WITHOUT regularization
        time_range = np.linspace(0,T_max, MEMORY_STEPS) # times for plotting
        
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))
        # Plot real part
        axs[0].plot(time_range, BCF.bcf(time_range).real, '-', label='Regularized')
        axs[0].plot(time_range, BCF_NOREG.bcf(time_range).real, '--', label='No Regularization')
        axs[0].set_xlabel('Time')
        axs[0].set_ylabel('Real BCF')
        axs[0].legend()

        # Plot imaginary part
        axs[1].plot(time_range, BCF.bcf(time_range).imag, '-', label='Regularized')
        axs[1].plot(time_range, BCF_NOREG.bcf(time_range).imag, '--', label='No Regularization')
        axs[1].set_xlabel('Time')
        axs[1].set_ylabel('Imaginary BCF')
        axs[1].legend()

        plt.tight_layout()
        plt.show()

    # TODO: figure out what this does
    itebd.gen_bath_from_BCF(BCF.bcf)
    
    itebd.compute_f(rtol=precision)
    
    # return the influence functional
    return itebd

#################### HAMILTONIAN FOR LANDAUER ERASURE ####################
# Time-dependent STA ramp function
def sta_ramp(t, tp, ramp_duration=1.0):
    """Smooth ramp function for STA activation."""
    if t < 0 or t > tp:
        return 0.0
    elif t < ramp_duration:
        # Ramp up: use raised cosine from 0 to 1
        return 0.5 * (1 - np.cos(np.pi * t / ramp_duration))
    elif t > tp - ramp_duration:
        # Ramp down: use raised cosine from 1 to 0
        return 0.5 * (1 - np.cos(np.pi * (tp - t) / ramp_duration))
    else:
        # Fully on in the middle
        return 1.0
    
    
# def H_ERASURE(t, tp, eps0, eps_max, n, shift=True, sta=False):
#     """
    
#     Calculate the Hamiltonian for Landauer erasure with a specified number of oscillations.
    
#     Parameters:
#         t (float): Current time.
#         tp (float): Protocolt time.
#         eps0 (float): Initial energy splitting.
#         epstau (float): Maximum energy splitting.
#         n (int): Number of oscillations (inverse) (1 for cyclic, 2 for Miller (2020) protocol).
#         shift (bool, optional): If True, apply a shift to the Hamiltonian.
#         sta (bool, optional): If True, apply the shortcut to adiavaticity (STA) term - so long as t>0. Default is False.
    
#    Returns:
#         h (numpy.ndarray): The Hamiltonian matrix.
   
#     """

#     # Calculate the mixing angle theta
#     theta = np.pi * ((t / tp) - 1)

#     # Calculate the time-dependent energy splitting
#     eps_of_t = eps0 + (eps_max - eps0) * (np.sin(np.pi * t / (n * tp)))**2

#     # Calculate the time dependent coefficients for sz and sx operators
#     fz = 0.5 * eps_of_t * np.cos(theta)
#     fx = 0.5 * eps_of_t * np.sin(theta)

#     # define the Hamiltonian
#     h  = fz * sz() + fx * sx() 
    
#     # if the shift is applied, add the shift term to the Hamiltonian
#     if shift:
#         # Calculate the shift required for ground state energy = 0 at all times
#         eta = np.sqrt(fz**2 + fx**2)
#         # add the shift to the hamiltonian with the identity operator
#         h += eta * (sz())**2

#     # if the shortcut to adiabaticity is applied, add the STA term
#     if sta:
#         if t>0 and t<tp: # apply it instantaneously at t=0 and turn it off at t=tp
#         # ramp = sta_ramp(t, tp) # ramp function for STA
#             ramp = 1 
#             h += ramp*(np.pi / (2 * tp)) * sy() 

#     # return the Hamiltonian
#     return h

def H_ERASURE(t, tp, eps0, eps_max, n, shift=True, sta=False, ramp_power=5):
    """
    Hamiltonian for smooth erasure protocol with tunable ramp shape.

    Parameters:
    - t: float, time
    - tp: float, total protocol duration
    - eps0: float, initial/final energy
    - eps_max: float, maximum energy in the middle
    - n: int, number of sin² oscillations in base modulation
    - shift: bool, add ground state energy shift
    - sta: bool, apply STA term (disabled by default)
    - ramp_power: float, power to raise the sin function for smoother turn-on/off
                  Higher values = flatter edges and steeper center
    """

    # Tunable ramp function: sin^ramp_power(pi * t / tp)
    ramp = np.sin(np.pi * t / tp)**ramp_power
    # Base energy modulation
    sin2_term = (np.sin(np.pi * t / (n * tp)))**2
    eps_of_t = eps0 + (eps_max - eps0) * sin2_term * ramp

    # Mixing angle for control direction
    theta = np.pi * ((t / tp) - 1)

    # Pauli terms
    fz = 0.5 * eps_of_t * np.cos(theta)
    fx = 0.5 * eps_of_t * np.sin(theta)

    # Hamiltonian construction
    h = fz * sz() + fx * sx()

    # Optional: zero ground state energy shift
    if shift:
        eta = np.sqrt(fz**2 + fx**2)
        h += eta * (sz()**2)

    # if the shortcut to adiabaticity is applied, add the STA term
    if sta:
        if t>0 and t<tp: # apply it instantaneously at t=0 and turn it off at t=tp
        # ramp = sta_ramp(t, tp) # ramp function for STA
            ramp = 1 
            h += ramp*(np.pi / (2 * tp)) * sy() 

    return h


#################### PULSES FOR WORK COUNTING ####################

def H_fw(j, s, m, f, step_size, tp, eps0, epstau, n, shift, sta):
    """
    Forward protocol Hamiltonian at discrete step j for work counting.

    Parameters:
        j (int): Current step index.
        s, m, f (int): Protocol step counts: equilibration, counting, and propagation steps.
        step_size (float): Time step size.
        tp (float): Total protocol time.
        eps0 (float): Initial energy splitting.
        epstau (float): Maximum energy splitting.
        n (int): n=1 for cyclic process, n=2 for half of the cycle.
        shift (bool): Whether to apply energy shift.
        sta (bool): Whether to apply shortcut to adiabaticity.

    Returns:
        numpy.ndarray: Hamiltonian matrix at step j during forward protocol.
    """

    # depending on the step j along the generalised time axis, use a different Hamiltonian
    if j <= s + m:
        t = 0
    else:
        t = (j - s - m) * step_size
    return H_ERASURE(t, tp, eps0, epstau, n, shift, sta)

def H_bw(j, s, m, f, step_size,tp, eps0, epstau, n, shift, sta):
    """
    Backward protocol Hamiltonian at discrete step j for work counting.

    Parameters:
        j (int): Current step index.
        s, m, f (int): Protocol step counts: equilibration, counting, and propagation steps.
        step_size (float): Time step size.
        tp (float): Total protocol time.
        eps0 (float): Initial energy splitting.
        epstau (float): Maximum energy splitting.
        n (int): n=1 for cyclic process, n=2 for half of the cycle.
        shift (bool): Whether to apply energy shift.
        sta (bool): Whether to apply shortcut to adiabaticity.

    Returns:
        numpy.ndarray: Hamiltonian matrix at step j during forward protocol.
    """
    # depending on the step j along the generalised time axis, use a different Hamiltonian
    if j <= s:
        t = 0
    elif j>s and j<=s+f:
        t = (j-s)  * step_size
    else:
        t = f * step_size
    
    return H_ERASURE(t, tp, eps0, epstau, n, shift, sta)

def Lio(j, s, m, f, step_size, tp, eps0, epstau, n, shift, sta):
    """
    Compute the Liouvillian superoperator at step j based on forward and backward Hamiltonians.

    Parameters:
        j, s, m, f, step_size, tp, eps0, epstau, n, shift, sta: Same as for H_fw and H_bw.

    Returns:
        numpy.ndarray: The Liouvillian operator matrix at step j.
    """
    
    # define the forward and backward Hamiltonians
    hfw = H_fw(j, s, m, f, step_size,tp, eps0, epstau, n, shift, sta)
    hbw = H_bw(j, s, m, f, step_size, tp, eps0, epstau, n, shift, sta)
    
    # define the Liouvillian superoperator
    liouvillian = -1j * (spre(hfw) - spost(hbw))

    return liouvillian

def expo_prop(j, s, m, f, step_size,tp, eps0, epstau, n, shift, sta):
    """
    Calculate the exponential propagator for the system at step j.

    Parameters:
        j, s, m, f, step_size, tp, eps0, epstau, n, shift, sta: Same as for Lio.

    Returns:
        numpy.ndarray: Matrix exponential of the Liouvillian operator scaled by half the step size.
    """
    return expm(Lio(j, s, m, f, step_size,tp, eps0, epstau, n, shift, sta) * step_size/2)

def U_tens(j, s, m, f, step_size, tp, eps0, epstau, n, shift, sta):
    """

    Parameters:
        j, s, m, f, step_size, tp, eps0, epstau, n, shift, sta: Same as for expo_prop.

    Returns:
        numpy.ndarray: A 3D tensor representing the propagator at step j.
    """
    liu_s = expo_prop(j, s, m, f, step_size, tp, eps0, epstau, n, shift, sta)
    u = np.einsum('ab,bc->abc', liu_s.T, liu_s.T)
    return u

#################### PROPAGATION OF GENERALISED DENSITY OPERATOR ####################

def process_iteration(itebd, s, f, m, rho0, STEP_SIZE, tp, eps0, epstau, n, shift, sta):
    """
    Perform a single iteration of system evolution using ITEBD with given propagators and initial state.

    Parameters:
        itebd (object): Influence functional object.
        s, m, f (int): Protocol step counts: equilibration, counting, and propagation steps.
        rho0 (numpy.ndarray): Initial density matrix.
        STEP_SIZE (float): Time step size.
        tp (float): Total protocol time.
        eps0 (float): Initial energy splitting.
        epstau (float): Maximum energy splitting.
        n (int): n=1 for cyclic process, n=2 for half of the cycle.
        shift (bool): Whether to apply energy shift.
        sta (bool): Whether to apply shortcut to adiabaticity.

    Returns:
    complex: Trace of the final density matrix (characteristic function) after evolution.
    """
    # create the list of propagators for the system stored contiguously in memory
    prop_list = np.ascontiguousarray([U_tens(j, s, m, f, STEP_SIZE, tp, eps0, epstau, n, shift, sta) 
                          for j in range(0, s + m + f + 1)])
    
    # extract the generalised density operator at the final time
    generalised_den_op = itebd.extract_dynamics(rho0, prop_list, only_final=True) # extract final state only

    # calculate the characteristic function as the trace of the final generalised density operator
    work_cf = np.trace(generalised_den_op.reshape(2, 2))

    # return the characteristic function
    return work_cf


def system_dynamics(itebd, s, f, m, rho0, STEP_SIZE,  tp, eps0, epstau, n , shift, sta):
    """
    Calculate the system states during the protocol.

    Parameters:
        itebd (object): Influence functional object.
        s, m, f (int): Protocol step counts: equilibration, counting, and propagation steps.
        rho0 (numpy.ndarray): Initial density matrix.
        STEP_SIZE (float): Time step size.
        tp (float): Total protocol time.
        eps0 (float): Initial energy splitting.
        epstau (float): Maximum energy splitting.
        n (int): n=1 for cyclic process, n=2 for half of the cycle.
        shift (bool): Whether to apply energy shift.
        sta (bool): Whether to apply shortcut to adiabaticity.

    Returns:
        states: A list of system states.
    """

    # create the list of propagators for the system, stored contiguously in memory
    prop_list = np.ascontiguousarray([U_tens(j, s, m, f, STEP_SIZE, tp, eps0, epstau, n , shift, sta) 
                          for j in range(0, s + m + f + 1)])
  
    # extract the states at all times using only_final=False
    states = itebd.extract_dynamics(rho0, prop_list, only_final=False) # extract all states
    
    # uncomment the line below to reshape the states into 2x2 matrices
    # states_matrix = [state.reshape(2, 2) for state in states] # return 2x2 matrices representing rho(t)
    
    return states

#################### EXTRACT SYSTEM PROPERTIES ####################

def trace_mag_coh(states):
    """
    Calculate the trace, magnetization, and coherence for a series of density matrices.

    Parameters:
        states (numpy.ndarray): Array of density matrices (each 2x2).

    Returns:
        tuple: Three numpy arrays containing traces, magnetizations (⟨sz⟩), and coherences (⟨sx⟩) at each state.
    """

    # Initialize lists to store trace, magnetization, and coherence
    tr = []
    mag = []
    coh = []

    for index in range(states.shape[0]): #loop through each state
        state = states[index] # 
        tr.append(np.trace(state)) # trace
        mag.append(np.trace(state.dot(sz())).real) # magnitization
        coh.append(np.trace(state.dot(sx())).real) # coherence

    # Convert lists to numpy arrays and return
    return np.array(tr), np.array(mag), np.array(coh)

def system_energy(states, tlist,  tp, eps0, epstau, n , shift, sta):
    """
    Compute the expectation value of the system energy for each state at given times.

    Parameters:
        states (numpy.ndarray): Array of density matrices.
        tlist (list or numpy.ndarray): List of time points corresponding to states.
        tp, eps0, epstau, n, shift, sta: Parameters for the Hamiltonian.

    Returns:
        numpy.ndarray: Energy expectation values at each time point.
    """

    # Ensure states and tlist have the same length
    if states.shape[0] != len(tlist):
        raise ValueError("The number of states and time points (tlist) must be the same.")
    
    # Initialize list to store energy values
    energy = [] 

    for index in range(states.shape[0]): # loop through number of states (same as number of timesteps)
        
        state = states[index] # pick out state 
        time = tlist[index] # pick out time

        # Calculate the energy expectation value at this time, taking the trace of the state with the Hamiltonian
        energy.append(np.trace(state.dot(H_ERASURE(time, tp, eps0, epstau, n , shift, sta))).real)
    
    # Convert the list of energy values to a numpy array and return
    return np.array(energy)

def overlap(states, tlist, tp, eps0, epstau, n , shift, sta):
    """
    Compute overlaps of system states with the instantaneous ground state and target ground state.

    Parameters:
        states (numpy.ndarray): Array of density matrices.
        tlist (list or numpy.ndarray): List of time points corresponding to states.
        tp, eps0, epstau, n, shift, sta: Parameters for the Hamiltonian.

    Returns:
        tuple: Two numpy arrays containing overlaps with the target ground state (at tp) and instantaneous ground states.
    """

    # Initialize lists to store overlaps
    overlap_target = []
    overlap_instant = []
    
    # find H(tp) ground state for overlap calculations
    H_tp =  H_ERASURE(tp, tp, eps0, epstau, n , shift, sta) # Hamiltonian at tp
    _, eigenvectors = np.linalg.eigh(H_tp) # eigenvectors of the Hamiltonian at tp
    ground = eigenvectors[0] # ground state at tp, 2x1 vector
    target = np.outer(ground, np.conj(ground)) # Outer product with its conjugate to make 2x2 matrix
    
    # loop through time
    for c, t in enumerate(tlist):   
        # pick out reduced system state at this time
        state = states[c] # 2x2 matrix  
        
        # calculate target overlap
        overlap_target.append(np.trace(target.dot(state)).real) # only store the real part
         
        # find instantaneous ground state at this time
        Ht =  H_ERASURE(t, tp, eps0, epstau, n, shift, sta) # hamiltonian at time t
        _, evec = np.linalg.eigh(Ht) # eigenvectors of the hamiltonian at time t
        gt = evec[0] # ground state at time t, 2x1 vector
        instant_ground = np.outer(gt, np.conj(gt)) # Outer product with its conjugate to make 2x2 matrix
        
        # calculate overlap with instantaneous ground state
        overlap_instant.append(np.trace(instant_ground.dot(state)).real) # only store the real part

    # Convert lists to numpy arrays and return
    return np.array(overlap_target), np.array(overlap_instant)

#################### OTHER USEFUL FUNCTIONS ####################

def get_deriv_coeffs(MAX_DIFF_ORDER, NTH_ORDER_DERIV):
    """
    Generate finite difference coefficients for derivatives up to fourth order with varying accuracies.

    Parameters:
        MAX_DIFF_ORDER (int): Maximum order of derivative accuracy (even numbers).
        NTH_ORDER_DERIV (int): Number of higher order derivatives to include beyond the base.

    Returns:
    tuple: Four numpy arrays containing finite difference coefficients for first, second, third, and fourth derivatives.
           Each array is shaped (MAX_DIFF_ORDER, MAX_DIFF_ORDER + NTH_ORDER_DERIV).
    """

    # Check if MAX_DIFF_ORDER is a valid even integer and NTH_ORDER_DERIV is non-negative
    if MAX_DIFF_ORDER < 2 or MAX_DIFF_ORDER % 2 != 0:
        raise ValueError("MAX_DIFF_ORDER must be an even integer ≥ 2.")
    if NTH_ORDER_DERIV < 0:
        raise ValueError("NTH_ORDER_DERIV must be non-negative.")

    # Initialize arrays for finite difference coefficients
    first_deriv_arr = np.zeros((MAX_DIFF_ORDER, MAX_DIFF_ORDER+NTH_ORDER_DERIV))
    second_deriv_arr = np.zeros((MAX_DIFF_ORDER, MAX_DIFF_ORDER+NTH_ORDER_DERIV))
    third_deriv_arr = np.zeros((MAX_DIFF_ORDER, MAX_DIFF_ORDER+NTH_ORDER_DERIV))
    fourth_deriv_arr = np.zeros((MAX_DIFF_ORDER, MAX_DIFF_ORDER+NTH_ORDER_DERIV))

    # these arrays have columns for each step along chi, and columns for each accuracy of the approximation
    # findiff only gives you even orders of accuracy, so every odd row is set to zero

    for i in np.arange(2, MAX_DIFF_ORDER+2, step=2): # the 2 is because we only get even orders of accuracy from findiff
        
        # first derivative coefficients
        coefs = findiff.coefficients(deriv=1, acc=i, symbolic=False)
        co = coefs['forward']['coefficients']
        offset = coefs['forward']['offsets']
        first_deriv_arr[i-1, offset] = co

        # second derivative coefficients
        coefs = findiff.coefficients(deriv=2, acc=i, symbolic=False)
        co = coefs['forward']['coefficients']
        offset = coefs['forward']['offsets']
        second_deriv_arr[i-1, offset] = co

        # third derivative coefficients
        coefs = findiff.coefficients(deriv=3, acc=i, symbolic=False)
        co = coefs['forward']['coefficients']
        offset = coefs['forward']['offsets']
        third_deriv_arr[i-1, offset] = co
        
        # fourth derivative coefficients
        coefs = findiff.coefficients(deriv=4, acc=i, symbolic=False)
        co = coefs['forward']['coefficients']
        offset = coefs['forward']['offsets']
        fourth_deriv_arr[i-1, offset] = co

    return  first_deriv_arr, second_deriv_arr, third_deriv_arr, fourth_deriv_arr

