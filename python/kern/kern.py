import sys
sys.path.append('/home/james/mlprojects/ndlutil/python')
import ndlutil
import numpy as np

class kern(ndlutil.parameterised):
	"""The base class for kernels. Should not be instantiated"""
	def __init__(self):
		self.shape = tuple()
		self.Nparam = 0
	def set_param(self):
		raise NotImplementedError
	def get_param(self):
		raise NotImplementedError
	def get_param_names(self):
		raise NotImplementedError
	def set_X(self):
		raise NotImplementedError
	def compute(self):
		raise NotImplementedError
	def compute_diag(self):
		raise NotImplementedError
	def gradients(self):
		raise NotImplementedError
	def gradients_X(self):
		raise NotImplementedError
	def __add__(self,other):
		if isinstance(other,compound):
			other.__add__(self)
		else:
			return compound([self,other])

class compound(kern):
	def __init__(self,kerns=[]):
		if isinstance(kerns,kern):
			kerns = [kern]

		self.kerns = []
		self.shape = kerns[0].shape
		self.Nparam = 0
		self.param_allocation = []
		[self.__add__(k) for k in kerns]

	def set_param(self,x):
		[k.set_param(x[i]) for k,i in zip(self.kerns, self.param_allocation)]
	def get_param(self):
		return np.hstack([k.get_param() for k in self.kerns])
	def get_param_names(self):
		return sum([['cmpnd_%i_%s'%(i,n) for n in k.get_param_names()] for i,k in enumerate(self.kerns)],[])
	def set_X(self,X):
		[k.set_X(X) for k in self.kerns]
	def compute(self):
		return np.sum([k.compute() for k in self.kerns],0)
	def diag(self):
		return np.sum([k.diag() for k in self.kerns],0)
	def gradients(self):
		return sum([k.gradients() for k in self.kerns],[])
	def __add__(self,other):
		if isinstance(other,compound):
			assert self.shape==other.shape
			self.kerns.extend(other.kerns)
			self.param_allocation.extend([i+self.Nparam for i in other.param_allocation])
			self.Nparam += other.Nparam
		elif isinstance(other,kern):
			self.kerns.append(other)
			assert self.shape==other.shape
			self.param_allocation.append(np.arange(self.Nparam,self.Nparam+other.Nparam))
			self.Nparam += other.Nparam

		else:
			raise AttributeError, "%s is not a proper kernel"%str(k)

def hierarchical(connections,kerns):
	"""
	Arguments
	---------
	connections : a numpy boolean array with as many rows as X.
	kerns : a list of kernels  withthe same domain X.

	Notes
	-----
	the n,k th element of connections is True if the nth input is affected by the kth kernel

	Example
	-------
	con = np.hstack((np.ones((90,1)),np.kron(np.eye(3),np.ones((30,1)))))
	X = np.linspace(0,10,90)[:,None]
	kerns = [kern.rbf(X) for i in range(4)]
	k = kern.hierarchical(con,kerns)

	"""
	for i,k in enumerate(kerns):
		assert k.shape==kerns[0].shape
		k.set_mask(np.nonzero(connections[:,i])[0])
	return compound(kerns)


