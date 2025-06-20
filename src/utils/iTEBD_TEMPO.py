# Standard library imports
from typing import Callable, Optional

# Third-party imports
import numpy as np
from scipy.linalg import expm, norm, svd
from scipy.sparse.linalg import eigs
from tqdm import tqdm
from ncon import ncon  # Efficient tensor contraction

# Local imports
from utils.operators import spre, spost, sprepost, lind
from utils.bath_functions import BathCorrelation


def tanh_cutoff(b, t, x):
    return (t - b) / 2 + b -  np.tanh(x - 2.5) * (t - b) / 2

def svd_filter(b, t, endtime, nc):
    y = [tanh_cutoff(t, b, i) for i in np.linspace(0, endtime, nc)]
    vals = 10**(-np.array(y[::-1])) #reverse array
    return vals

def iTEBD_apply_gate(gate: np.ndarray, 
                     A: np.ndarray, 
                     sAB: np.ndarray, 
                     B: np.ndarray, 
                     sBA: np.ndarray, 
                     rank: int, 
                     rtol: float, 
                     ctol: Optional[float] = 1e-13):
    """
    single iTEBD step, scheme adapted from https://www.tensors.net/mps
    :param gate: TEBD gate for A-B link
    :param A: A tensor (left)
    :param sAB: weight for A-B link
    :param B: B tensor (right)
    :param sBA: weight for B-A link
    :param rank: maximum rank in svd compression
    :param rtol: relative error for svd compression
    :param ctol: cutoff for weights sBA which need to be inverted
    :return: new tensors and weights A, sAB, B, sBA
    """

    # renormalize weights
    sAB = sAB * norm(sBA)
    sBA = sBA / norm(sBA)

    # ensure weights are above tolerance (needed for inversion)
    #sBA[np.abs(sBA) < ctol] = ctol
    

    # MPS - gate contraction
    d1 = gate.shape[1]
    d2 = gate.shape[-1]
    rank_BA = sBA.shape[0]
    contractions = [[-1, 1], [1, 5, 2], [2, 4], [4, 6, 3], [3, -4], [5, -2, 6, -3]]
    tens = ncon([np.diag(sBA), A, np.diag(sAB), B, np.diag(sBA), gate], contractions)
    tens = tens.reshape([d1 * rank_BA, d2 * rank_BA])
    u, s_vals, v =svd(tens, full_matrices=False, lapack_driver='gesvd')

    # truncate singular values
    if rtol is None:
        chi = min(rank, len(s_vals))
    else:
        
        # s_vals_sum = np.cumsum(s_vals) / np.sum(s_vals)
        # rank_rtol = np.searchsorted(s_vals_sum, 1 - rtol) + 1
        # chi = min(rank, len(s_vals), rank_rtol)
        chi = (s_vals>rtol*np.max(s_vals)).sum()
    

    u = u[:, :chi].reshape(sBA.shape[0], d1 * chi)
    v = v[:chi, :].reshape(chi * d2, rank_BA)

    # factor out sAB weights from A and B
    A = (np.diag(1 / sBA) @ u).reshape(sBA.shape[0], d1, chi)
    B = (v @ np.diag(1 / sBA)).reshape(chi, d2, rank_BA)

    # new weights
    sAB = s_vals[:chi]

    return A, sAB, B, sBA




class iTEBD_TEMPO():
    """ A class to compute and approximate the influence functional using iTEBD-TEMPO and compute dynamics. """

    def __init__(self, 
                s_vals: np.ndarray, 
                STEP_SIZE: float, 
                MAX_STEPS: int):
        """
        :param s_vals: Real eigenvalues of the system-bath coupling operator.
        :param STEP_SIZE: time step for Trotter splitting.
        :param MAX_STEPS: Memory cutoff. Should be chosen large enough.
        """
        print("Include if statement...")
        self.bath = None
        self.n_c = MAX_STEPS

        self.s_vals = s_vals
        self.s_dim = self.s_vals.size
        self.nu_dim = self.s_vals.size ** 2 + 1
        
        self.STEP_SIZE = STEP_SIZE

        #define the commutator and anticommutator in Lioville space
        self.Ocomm = spre(np.diag(s_vals))-spost(np.diag(s_vals))
        self.Oanti = spre(np.diag(s_vals)) + spost(np.diag(s_vals))

        # extract the diagonal elements of the commutator and anticommutator
        self.s_sum = self.Oanti.diagonal()
        self.s_diff = self.Ocomm.diagonal()
        
        #add an extra dimension for the auxiliary space:
        self.s_sum = np.pad(self.s_sum, [(0, 1)])
        self.s_diff = np.pad(self.s_diff, [(0, 1)])
        
        self.kron_delta = np.identity(self.nu_dim)
        self.f = None
        print('ctol turned off...')


    def gen_bath_from_spectral_density(self, 
                Spectral_Function: Callable[[float], complex], 
                beta: float, 
                cutoff: Optional[float] = np.inf):
        """
        Generate the bath correlation function.

        :param Spectral_Function: Spectral function of the bath.
        :param beta: Inverse temperature of the bath.
        :param cutoff: Cutoff frequency for the bath.
        :return: Bath correlation function.
        """
        self.bath = BathCorrelation()
        self.bath.construct_eta_list( 
                           self.STEP_SIZE, 
                           self.n_c,
                           Spectral_Function,
                           beta, 
                           cutoff)
        

    def gen_bath_from_BCF(self, 
                BCF: Callable[[float], complex]):
        """
        Generate the bath correlation function.

        :param Spectral_Function: Spectral function of the bath.
        :param beta: Inverse temperature of the bath.
        :param cutoff: Cutoff frequency for the bath.
        :return: Bath correlation function.
        """
        self.bath = BathCorrelation()
        self.bath.set_bath_correlation(BCF)
        self.bath.compute_eta_from_BCF( 
                        self.n_c,
                        self.STEP_SIZE)

   
    def gen_gate(self, k: int)-> np.ndarray:
        i_tens = np.exp(-self.bath.eta[k].real * np.outer(self.s_diff, self.s_diff) 
                        - 1j * self.bath.eta[k].imag * np.outer(self.s_sum, self.s_diff))

        if k == 0:
            gate = np.einsum('a,ij,jb,j->jabi', 
                            np.ones((1)), 
                            self.kron_delta, 
                            self.kron_delta, 
                            np.diagonal(i_tens))
        else:
            gate = np.einsum('ij,ab,aj->jabi',  
                            self.kron_delta, 
                            self.kron_delta, 
                            i_tens)
        return gate

    def compute_f(self, rtol: float, rank: Optional[int] = np.inf,
                  ctol: Optional[float] = 1e-13):
        """
        Compute the infinite influence functional tensor f using iTEBD.

        :param rtol: Relative tolerance for svd compression.
        :param rank: Maximum allowed rank (bond dimension).
        """

        A = np.ones((1, self.nu_dim, 1))
        B = np.ones((1, self.nu_dim, 1))
        sAB = np.ones((1))
        sBA = np.ones((1))
        self.rank = []
        self.max_svd = []
        self.min_svd = []
        rank_is_one = True

        # rtols = svd_filter(5, 7.5, 6, self.n_c)

        for k in tqdm(range(1, self.n_c + 1), desc='building influence functional'):
            
            gate = self.gen_gate(self.n_c - k)

            if k % 2 == 0:
                # B, sBA, A, sAB = iTEBD_apply_gate(gate, B, sBA, A, sAB, rank, rtol=rtols[k-1],ctol=ctol)
                B, sBA, A, sAB = iTEBD_apply_gate(gate, B, sBA, A, sAB, rank, rtol=rtol,ctol=ctol)
                
                self.max_svd.append(max(sBA))
                self.min_svd.append(min(sBA))

            else:
                # A, sAB, B, sBA = iTEBD_apply_gate(gate, A, sAB, B, sBA, rank, rtol=rtols[k-1], ctol=ctol)
                A, sAB, B, sBA = iTEBD_apply_gate(gate, A, sAB, B, sBA, rank, rtol=rtol, ctol=ctol)
                
                self.max_svd.append(max(sAB))
                self.min_svd.append(min(sAB))
            self.rank.append(max(sAB.size, sBA.size))
            # print(self.rank[-1])
            if rank_is_one:
                if np.all([sAB.shape[0] == 1, sAB.shape[-1] == 1, sBA.shape[0] == 1, sBA.shape[-1] == 1]):
                    # reset to initial mps if rank is still one
                    sAB = np.ones((1))
                    sBA = np.ones((1))
                    A = np.ones((1, self.nu_dim, 1))
                    B = np.ones((1, self.nu_dim, 1))
                else:
                    rank_is_one = False
                    self.n_c_eff = self.n_c - k + 1
                    if k == 1:
                        print("""Warning: the memory cutoff n_c may be too small 
                              for the given rtol value. The algorithm may become unstable and inaccurate. 
                              It is recommended to increase n_c until this message does no longer appear.""")
        contractions=[[-1, 1], [1, -2, 2], [2, 3], [3, -3, -4]]
        
        
        # self.f = np.squeeze(ncon([np.diag(sAB), B, np.diag(sBA), A], contractions))

        AA = ncon([np.diag(sAB), B, np.diag(sBA), A], contractions) 
        self.f = AA[: :, 0, :] # allows alpha=0

        # compute the endcaps of the influence functional
        w, v_r = eigs(self.f[:, -1, :], 1, which='LR')
        w, v_l = eigs(self.f[:, -1, :].T, 1, which='LR')
        self.v_r = v_r[:, 0]
        self.v_l = v_l[:, 0] / (v_l[:, 0] @ v_r[:, 0])

        return

    def extract_single_time_maps(self, 
            h_s: np.ndarray, 
            n: int,
            dissipator=[]) -> np.ndarray:
        """
        Compute the maps for n time steps.

        :param h_s: System Hamiltonian in the eigenbasis of the coupling operator.
        :param n: Number of time-steps for the propagation.
        :return: Time evolution of density matrix.
        """

        assert self.f is not None, "the influence functional has not yet been computed, run self.compute_f first"
        lio = -1j *(spre(h_s) - spost(h_s))
        if len(dissipator) != 0:
            for rate, op in dissipator:

                lio += 0.5 * rate * lind(op)
                
        liu_s = expm(lio * self.STEP_SIZE / 2)
        u = np.einsum('ab,bc->abc', liu_s.T, liu_s.T)

        # maps_t = np.empty((n + 1, self.nu_dim-1, self.nu_dim-1), dtype=np.complex128)
        self.single_time_maps =  np.empty((n , self.nu_dim-1, self.nu_dim-1), dtype=np.complex128)
        

        evol_tens = ncon([self.f[:, :-1, :], u], [[-1, 2, -3], [-2, 2, -4]])
        # state = ncon([self.v_l, rho_0.flatten()], [[-2], [-3]])

        m = ncon([self.v_l, evol_tens], [[1], [1, -1, -2, -3]])
        # print(m.shape)
        self.single_time_maps[0] = ncon([m, self.v_r], [[-1, 1, -2], [1]])
        for i in range(1,n):
            m = ncon([m, evol_tens], [[-1,1,2],[1,2,-2,-3]])
            self.single_time_maps[i] = ncon([m, self.v_r], [[-1, 1, -2], [1]])
