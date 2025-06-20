# Standard library imports
from typing import Callable, Optional  

# Third-party imports
import numpy as np  
from scipy.integrate import quad, dblquad  
from scipy.interpolate import interp1d  
from tqdm import tqdm  
import matplotlib.pyplot as plt  

class BathCorrelation():
    """Class to represent and manipulate the bath correlation function."""

    def __init__(self):
        """Initializes the BathCorrelation instance."""

        self.bcf = None
        self.eta = None

    def set_bath_correlation(self, bcf: Callable[[float], complex]):
        """Sets the bath correlation function directly.

        Args:
            bcf: A callable that takes a time (float) and returns a complex value.
        """
        self.bcf = bcf
        self.time_cutoff = np.inf

    def construct_bath_correlation(self, 
                                   Spectral_Function: Callable[[float], float],
                                   beta: float,
                                   time_cutoff: float,
                                   time_steps: int,
                                   cutoff: Optional[float] = np.inf,
                                regularize: Optional[list] = [False, None]):
        """Constructs the bath correlation function from a spectral function.

        Args:
            spectral_function: The spectral function of the bath.
            beta: The inverse temperature.
            time_cutoff: The cutoff time for the correlation function.
            time_steps: Number of time steps for discretization.
            cutoff: Cutoff frequency for the spectral function.
            regularize: Optional list where first element is a boolean to apply tapering,
                and second element is a tapering function.
        
        Sets:
            self.bcf: An interpolating function valid on [0, time_cutoff].
        """

        re_int = lambda t, omega:  Spectral_Function(omega) * np.cos(omega * t) /np.tanh(0.5 * beta * omega)
        im_int = lambda t, omega: Spectral_Function(omega) * np.sin(omega * t)
        bcf_re_data = []
        bcf_im_data = []

        trange = np.linspace(0, time_cutoff, time_steps)
        
        for t in tqdm(trange, desc="Constructing bath correlation function."):
            if regularize[0]:

                taper_func = regularize[1](t)
            else: 
                taper_func = 1

            re,_ = quad(lambda omega: re_int(t, omega), 0, cutoff,limit=100000)
            im,_ = quad(lambda omega: im_int(t, omega), 0, cutoff,limit=100000)
            bcf_re_data.append(re*taper_func)
            bcf_im_data.append(im*taper_func)
        
        interp_re = interp1d(trange, bcf_re_data)
        interp_im = interp1d(trange, bcf_im_data)
        
        
        self.bcf = lambda t: interp_re(t) - 1j * interp_im(t)
        self.time_cutoff = time_cutoff

    def numerical_eta(self,t: float, 
                      Spectral_Function: Callable[[float], float],
                        beta: float, 
                        cutoff: Optional[float] = np.inf) -> complex:
        """Evaluates the memory kernel at time t via numerical integration.

        The time integral is performed analytically, returning the memory kernel.

        Args:
            t: Time at which to evaluate the memory kernel.
            spectral_function: The spectral function of the bath.
            beta: The inverse temperature.
            cutoff: Cutoff frequency for the spectral function.

        Returns:
            Complex value of the memory kernel at time t.
        """
        re_func = lambda w: (Spectral_Function(w)/w**2) * (1-np.cos(w * t)) / np.tanh(0.5 * w * beta)
        #re_data = * coth_list
        im_func = lambda w: (Spectral_Function(w)/w**2) * (np.sin(w * t) - w * t)

        re,_ = quad(re_func, 0, cutoff, limit=1000)
        im,_ = quad(im_func, 0, cutoff, limit=1000)
        return re + 1j * im
    
    def construct_eta_list(self, 
                           step_size: float, 
                           number_of_steps: int,
                           Spectral_Function: Callable[[float], float],
                           beta: float, 
                           cutoff: Optional[float] = np.inf):
        """Constructs the memory kernel (eta) over a range of times.

        Args:
            step_size: Time step size.
            number_of_steps: Number of time steps.
            spectral_function: The spectral function of the bath.
            beta: The inverse temperature.
            cutoff: Cutoff frequency for the spectral function.
        
        Sets:
            self.eta: Complex numpy array of memory kernel values.
        """
        eta_list = []
        for n in range(number_of_steps+1):
            eta_list.append(self.numerical_eta(n * step_size, 
                                               Spectral_Function, 
                                               beta, cutoff))
        print(eta_list[0], eta_list[1])
        self.eta = np.empty(number_of_steps, dtype=np.complex128)
        for k in range(number_of_steps):
            if k==0:
                self.eta[k] = eta_list[k+1]
            else:
                self.eta[k] = eta_list[k + 1] - 2 * eta_list[k] + eta_list[k-1]

    def check_bath_correlation(self,t_end: float, time_steps: int):
        """Plots the real and imaginary parts of the bath correlation function.

        Args:
            t_end: End time for plotting.
            time_steps: Number of time steps in the plot.
        
        Raises:
            AssertionError: If bath correlation function is not set.
        """
        assert self.bcf is not None, "Bath correlation function not set."
        trange = np.linspace(0, t_end, time_steps)

        plt.plot(trange, np.real(self.bcf(trange)), '.-',label="Re")
        plt.plot(trange, np.imag(self.bcf(trange)),  '.-',label="Re")
        plt.show()
        plt.close()

    def compute_eta_from_BCF(self, NO_TIME_STEPS: int, TIME_STEP: float) -> np.ndarray:
        """Computes the discretized memory kernel (eta) from the bath correlation function.

        Args:
            no_time_steps: Number of discrete time steps.
            time_step: Size of each time step.

        Returns:
            Complex numpy array representing the discretized memory kernel.
        
        Raises:
            AssertionError: If bath correlation function is not set or time cutoff is too small.
        """
        assert self.bcf is not None, "Bath correlation function not set."
        assert self.time_cutoff > NO_TIME_STEPS * TIME_STEP, "Time cutoff is too small"

        self.eta = np.zeros(NO_TIME_STEPS, dtype=np.complex128)
        self.eta[0] += dblquad(lambda s, t: np.real(self.bcf(t - s)), 0, TIME_STEP, lambda t: 0, lambda t: t)[0]
        self.eta[0] += dblquad(lambda s, t: np.imag(self.bcf(t - s)), 0, TIME_STEP, lambda t: 0, lambda t: t)[0] * 1j

        for k in tqdm(range(1, NO_TIME_STEPS), desc="Computing eta function."):
            self.eta[k] += dblquad(lambda s, t: np.real(self.bcf(t - s)), k * TIME_STEP, (k + 1) * TIME_STEP, 0, TIME_STEP)[0]
            self.eta[k] += dblquad(lambda s, t: np.imag(self.bcf(t - s)), k * TIME_STEP, (k + 1) * TIME_STEP, 0, TIME_STEP)[0] * 1j
