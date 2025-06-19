# Third-party imports
import numpy as np 
from numpy import tensordot  

# Local project imports
from utils.iTEBD_TEMPO import iTEBD_TEMPO  

class time_dependent_TEMPO(iTEBD_TEMPO):

    def __init__(self, 
                 s_vals: np.ndarray, 
                 STEP_SIZE: float, 
                 MAX_STEPS: int):
        
        """
        :param s_vals: The eigenvalues of the coupling operator.
        :param STEP_SIZE: The time step size.
        :param MAX_STEPS: The number of timesteps.
        """
        #instantiating the iTEBD parent class         
        super().__init__(s_vals, STEP_SIZE, MAX_STEPS)

        # set the system propagator to None
        self.sys_prop = None

    def extract_dynamics(self,
                        Rho0: np.ndarray,
                        propagators: np.ndarray,
                        only_final = True) -> np.ndarray:
        """
        Extract the system dynamics using the computed influence functional.

        This method evolves the initial density matrix (Rho0) using the precomputed
        influence functional and the provided system propagators. It can either
        return the full time evolution of the system or just the final state.

        i.e. it returns rho(chi, t), from which we can calcualte the CF or moments

        Parameters:
        -----------
        Rho0 : np.ndarray
            The initial density matrix of the system.
        propagators : np.ndarray
            Array of system propagators for each time step.
        only_final : bool, optional
            If True, only the final state is returned. If False, the full time
            evolution is returned. Default is True.

        Returns:
        --------
        np.ndarray
            If only_final is True, returns the final density matrix.
            If only_final is False, returns an array of density matrices for each time step.

        Raises:
        -------
        AssertionError
            If the influence functional has not been computed before calling this method.
        """
        assert self.f is not None, "the influence functional has not yet been computed, run self.compute_f first"
        
        # Contract the influence functional with the first propagator
        # along axis=1 of each the IF and the first propagator
        evol_tens = tensordot(self.f[:, :-1, :], propagators[0], axes=[[1],[1]]) # remove last element of axis=1 for IF
        
        # Reorder the axes of the resulting tensor
        evol_tens = np.moveaxis(evol_tens, [0,1,2,3], [0,2,1,3])
     
        # Contract the left boundary vector with the evolution tensor
        m = tensordot(self.v_l, evol_tens, axes=[[0], [0]])
     
        # Contract the result with the flattened initial density matrix
        m = tensordot(m, Rho0.reshape(-1), axes=[[0],[0]])
      
        # Extract the raw influence functional (excluding the last time step)
        raw_f = self.f[:, :-1, :] # this is repeated above, can we save it?
        if only_final == False:
            # Initialize an array to store all states if full evolution is requested
            states = np.zeros((propagators.shape[0], Rho0.shape[0], Rho0.shape[1]), dtype = np.complex128)
            states[0] = Rho0  # Set the initial state
        
        for i in range(1, propagators.shape[0]):  # loop all timesteps except t=0
                    
            u = propagators[i]  # Get the current propagator
            
            # Contract m with the current propagator
            m = tensordot(m, u, axes=[[1],[0]])
            
            # Contract m with the raw influence functional
            m = tensordot(m,raw_f, axes=[[0,1],[0,1]]).T

            if only_final == False:
                # If full evolution is requested, store the current state
                states[i] = tensordot(m, self.v_r, axes=[[0],[0]]).reshape(Rho0.shape)
        
        if only_final == False:
            return states  # Return all states if full evolution was requested
        
        else:
            # Calculate and return only the final state
            state = tensordot(m, self.v_r, axes=[[0],[0]])
            return state
