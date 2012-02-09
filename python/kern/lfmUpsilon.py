import lfmUpsilonf2py as lfmulocal
import numpy as np

def compupvec(gamma, sig2, t1):
    """This function computes Upsilon when it is needed as a vector.
    It is called inside the lfm kernel. It does the same that
    lfmComputeUpsilonVector.m in the matlab version
    """
    return(lfmulocal.compupvec(gamma,sig2,t1,np.shape(t1)[0]))

def gradupvec(gamma, sig2, t1):
    """This function computes the gradient of Upsilon when it is needed as a vector.
    It is called inside the lfm kernel. It does the same that
    lfmGradientUpsilonVector.m in the matlab version.
    """
    return(lfmulocal.gradupvec(gamma,sig2,t1,np.shape(t1)[0]))

def compupmat(gamma, sig2, t1, t2):
    """This function computes Upsilon when it is needed as a matrix.
    It is called inside the lfm kernel. It does the same that
    lfmComputeUpsilonMatrix.m in the matlab version
    """
    return(lfmulocal.compupmat(gamma,sig2,t1,t2, np.shape(t1)[0], np.shape(t2)[0]))

def gradupmat(gamma, sig2, t1, t2):
    """This function computes the gradient of the Upsilon when it
    is needed as a matrix, It is called inside the lfm kernel. It does the same that
    lfmGradientUpsilonMatrix.m in the matlab version
    """
    return(lfmulocal.gradupmat(gamma,sig2,t1,t2, np.shape(t1)[0], np.shape(t2)[0]))

def gradsigupvec(gamma, sig2, t1):
    """This function computes the gradient of Upsilon wrt sigma when it is needed as a vector.
    It is called inside the lfm kernel. It does the same that
    lfmGradientSigmaUpsilonVector.m in the matlab version
    """
    return(lfmulocal.gradsigupvec(gamma,sig2,t1,np.shape(t1)[0]))

def gradsigupmat(gamma, sig2, t1, t2):
    """This function computes the gradient of Upsilon wrt sigma when it is needed as a matrix.
    It is called inside the lfm kernel. It does the same that
    lfmGradientSigmaUpsilonMatrix.m in the matlab version
    """
    return(lfmulocal.gradsigupmat(gamma,sig2,t1,t2, np.shape(t1)[0], np.shape(t2)[0]))

def h(gamma1, gamma2, sig2, t1, t2):
	"""
	Compute the function h
	"""
	return 	(compupmat(gamma1,sig2,t1,t2) + np.exp(-gamma2*t1)*compupmat(gamma,sig1,t1,0.*t2))/(gamma1 + gamma2)



