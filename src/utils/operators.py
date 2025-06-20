# Third-party imports
import numpy as np 
import scipy as sc 


def spre(op):
    """ Superoperator for the left action (pre) of an operator on a density matrix. """
    return np.kron(op, np.eye(op.shape[0]))


def spost(op):
    """ Superoperator for the right action (post) of an operator on a density matrix. """
    return np.kron(np.eye(op.shape[0]), op.T)


def sprepost(A, B):
    """ Superoperator for the left action of A and the right action of B on a density matrix. """
    #the Liouville space representation of A.R.B
    return np.kron(A, B.T)


def lind(op):
    """ Returns the Lindblad superoperator for a given operator `op`. """
    super_op = 2 * sprepost(op,op.T.conj())                        
    super_op -= spre(op.T.conj()@op)
    super_op -= spost(op.T.conj()@op)
    return super_op


def coth(x):
    """ Computes the hyperbolic cotangent of x, handling divergence. """
    try:
        val = np.cosh(x)/np.sinh(x)
    except RuntimeWarning:
        print('divergence...')
        val=0
    return val


def scipy_compress(A, prec):
    """
    This function finds the singular value decomposition for an arbitrary tensor
    A, and truncates up to some precision prec.
    """
    U, S, V= sc.linalg.svd(A, full_matrices=False, 
                           compute_uv=True, lapack_driver='gesvd') 
    
    chi = np.sum(S>prec*S[0])
    # print(len(S), chi)
    return chi, U[:,0:chi], S[0:chi], V[0:chi, :]


def sz():
    """ Returns the Pauli Z operator. """
    return np.array([[1,0],[0,-1]],dtype=np.complex128)

def sx():
    """ Returns the Pauli X operator. """
    return np.array([[0,1],[1,0]],dtype=np.complex128)

def sy():
    """ Returns the Pauli Y operator. """
    return np.array([[0,-1j],[1j,0]],dtype=np.complex128)



def calculate_indistinguishabilty(corr_func:np.ndarray,
                                  STEP_SIZE: float):
    """
    A function to calculate the indistinguishability of a two-time 
    correlation function.
    The input is a two-time correlation function with the first axis being the 
    tau axis and the second axis being the time axis. 
    returns the indistinguishability of the correlation function.  
    """    
    number_of_maps = corr_func.shape[1]
    
    #calculate the numerator:    
    num = (np.abs(corr_func)**2).sum() * STEP_SIZE**2

    #Calculate the denominator:
    #first pull out the single time correllation funciton:
    f = corr_func[0,:]

    #iteratively calculate the denominator:
    den = 0
    for n in range(number_of_maps):
        for j in range(n, number_of_maps):
            den += f[n] * f[j] * STEP_SIZE**2

    return num/den