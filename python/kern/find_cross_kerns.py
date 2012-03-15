
from lfm import *
from rbf import *

table = {
	lfm: {lfm:lfmXlfm, rbf:lfmXrbf},
	rbf: {lfm:rbfXlfm}
	}
#TODO: more kernels to follow. Sim, in particular
	
def find_cross_kerns(k1, k2):
	"""
	A helper function to find the appropriate cross-kernel object for a pair of kernels
	"""
	try:
		xk = table[k1][k2]
	except KeyError:
		raise AttributeError, "cross kernel not found"
	return xk(k1, k2)

