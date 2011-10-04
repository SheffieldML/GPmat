import numpy as np
from kern import kern

class stationary(kern):
	"""
	Base class for kernels which are a functions of the absolute difference r.
	See Rasmussen and Williams p 83.

	Notes
	-----
	To implement a new stationary covariance, simply define the normalised 
	function and its gradient wrt gamma. This class takes care of the rest!
	"""
	def __init__(self,X,alpha=1., gamma=1.):
		kern.__init__(self)
		self.masked=False
		self.alpha = alpha #variance
		self.gamma = gamma #parameter of function, to be defined by inherreted class ( usually a lengthscale)
		self.Nparam = 2
		self.set_X(X)

	def set_X(self,X):
		self.r = np.sum(np.abs(X[:,None,:]-X[None,:,:]),-1)
		self.X = X
		self.shape = self.r.shape
		if self.masked:
			self.r = self.r[self.mask]

	def set_mask(self,i):
		"""
		Apply a mask to this kernel. i is a list or array of integer indices
		"""
		self.masked=True
		self.mask = np.meshgrid(i,i)
		self.set_X(self.X)

	def compute(self,target=None):
		"""
		Arguments
		---------
		target : a np array to add the computation of this kernel to. If None, a new array is created
		"""
		if target is None:
			target = np.zeros(self.shape)
		if self.masked:
			target[self.mask] += self.alpha*self.function(self.r)
		else:
			target[:,:] += self.alpha*self.function(self.r)
		return target
		
	def diag_compute(self):
		if self.masked:
			ret = np.zeros(self.shape[0])
			ret[self.mask[0][0]] = self.alpha
			return ret
		else:
			return np.ones(self.shape[0])*self.alpha

	def cross_compute(self,X2):
		if self.masked:
			ret = np.zeros((self.shape[0],X2.shape[0]))
			r = np.sum(self.X[self.mask[0][0]][:,None,:]-X2[None,:,:],-1)
			ret[self.mask[0][0],:] = self.alpha*self.function(r)
			return ret
		else:
			r = np.square(np.sum(self.X[:,None,:]-X2[None,:,:],-1))
			return self.alpha*self.function(r)

	def function(self):
		raise NotImplementedError
	def function_gradients(self):
		raise NotImplementedError
	def gradients(self):
		return [self.function(), self.alpha*self.function_gradient()]
	def get_param(self):
		return np.array([self.alpha, self.gamma])
	def set_param(self,x):
		self.alpha,self.gamma = x
	def get_param_names(self):
		return ['alpha','gamma']

class rbf(stationary):
	"""
	Notes
	-----
	$$
	k(x,x') = \\alpha \\exp\{-\\gamma r^2\\}
	$$
	where $r = |x-x'|$
	"""
	def __init__(self,X,alpha=1., gamma=1.):
		stationary.__init__(self,X,alpha,gamma)
	def function(self,r):
		return np.exp(-self.gamma*np.square(r))
	def function_gradient(self,r):
		return -np.square(r)*self.function(r)

class Matern32(stationary):
	"""
	Notes
	-----
	$$
	k(x,x') = \\alpha (1+\\gamma r)\\exp\{-\\gamma r\\}
	$$
	where $r = |x-x'|$
	"""
	def __init__(self,X,alpha=1., gamma=1.):
		stationary.__init__(self,X,alpha,gamma)
	def function(self,r):
		return (1. + self.gamma*self.r)*np.exp(-self.gamma*self.r)
	def function_gradient(self):
		return -self.gamma*np.square(self.r)*np.exp(-self.gamma*self.r)

class Ornstein_Uhlenbeck(stationary):
	"""
	Notes
	-----
	$$
	k(x,x') = \\alpha \\exp\{-\\gamma r\\}
	$$
	where $r = |x-x'|$
	"""
	def __init__(self,X,alpha=1., gamma=1.):
		stationary.__init__(self,X,alpha,gamma)
	def function(self,r):
		return np.exp(-self.gamma*self.r)
	def function_gradient(self):
		return -self.gamma*self.r*np.exp(-self.gamma*self.r)











	
