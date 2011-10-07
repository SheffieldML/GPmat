import numpy as np
from kern import kern

class stationary(kern):
	"""
	Base class for kernels which are functions of the absolute difference r.
	See Rasmussen and Williams p 83.

	Notes
	-----
	To implement a new stationary covariance, simply define the normalised 
	function, its gradient wrt gamma and gradient wrt r. This class takes care of the rest!

	This also allows for efficient computation: r is stored when set_X is called.

	TODO
	----
	Perhaps some of the masking functionality could be moved to the kernel base class?
	"""
	def __init__(self,X,alpha=1., gamma=1.):
		kern.__init__(self)
		self.alpha = alpha #variance
		self.gamma = gamma #parameter of function, to be defined by inherreted class ( usually a lengthscale)
		self.Nparam = 2
		self.set_X(X)

	def set_X(self,X):
		"""Set self X  to the passed value and compute and store r"""
		self.r = np.sqrt(np.sum(np.square(X[:,None,:]-X[None,:,:]),-1))
		self.X = X
		self.shape = self.r.shape
		if self.masked:
			self.r = self.r[self.mask_grid]

	def compute(self,target=None):
		"""
		Arguments
		---------
		target : a np array to add the computation of this kernel to, for in-place computation. If None, a new array is created
		"""
		if target is None:
			target = np.zeros(self.shape)
		if self.masked:
			target[self.mask_grid] += self.alpha*self.function(self.r)
		else:
			target[:,:] += self.alpha*self.function(self.r)
		return target
		
	def diag_compute(self):
		"""Compute just the diagonal terms of the covariance"""
		if self.masked:
			ret = np.zeros(self.shape[0])
			ret[self.mask] = self.alpha
			return ret
		else:
			return np.ones(self.shape[0])*self.alpha

	def cross_compute(self,X2,target=None):
		"""Compute the covariance between self.X and the passed values X2. Mask aware"""
		if target is None:
			target = np.zeros((self.shape[0],X2.shape[0]))
		if self.masked:
			r = np.sqrt(np.sum(np.square(self.X[self.mask][:,None,:]-X2[None,:,:]),-1))
			target[self.mask,:] += self.alpha*self.function(r)
		else:
			r = np.sqrt(np.sum(np.square(self.X[:,None,:]-X2[None,:,:]),-1))
			target += self.alpha*self.function(r)
		return target

	def function(self):
		raise NotImplementedError
	def function_gradient(self):
		raise NotImplementedError
	def function_gradient_r_over_r(self):
		raise NotImplementedError
	def gradients(self):
		"""
		The gradient of the kernel wrt the variance and gamma.

		Returns
		-------
		A list of np arrays.
		"""
		if self.masked:
			ret = [np.zeros(self.shape), np.zeros(self.shape)]
			ret[0][self.mask_grid] = self.function(self.r)
			ret[1][self.mask_grid] = self.alpha*self.function_gradient(self.r)
		else:
			ret = [self.function(self.r), self.alpha*self.function_gradient(self.r)]
		return ret
	def gradients_X(self,target=None):
		"""
		The gradient of the kernel wrt the input points self.X

		In-place computation available by passing the target array. 

		Relies on the implementation of function_gradient_r_over_r:
		$$
		\\frac{\\partial f(r_{i,j}}{\\partial r_{i,j} \\times \\frac{1}{r_{i.j}}
		$$
		which arises from applying the chain rule to the derivative. Definin it this way helps prevent singlularity.

		Returns
		-------
		A NxNxD array. the i,j,kth element is the gradient of the [i,j]the element of the kernel wrt the k'th dimension of the i'th input.
		"""
		if target is None:
			target = np.zeros(self.shape+(self.X.shape[1],))
		if self.masked:
			common = np.zeros(self.shape,dtype=np.float64)
			common[self.mask_grid] = self.alpha*self.function_gradient_r_over_r(self.r)
			[np.add(target[:,i,:],common[i,:][:,None]*(self.X[i]-self.X),target[:,i,:]) for i in self.mask]
		else:
			common = self.alpha*self.function_gradient_r_over_r(self.r)
			[np.add(target[:,i,:],common[i,:][:,None]*(self.X[i]-self.X),target[:,i,:]) for i in range(self.shape[0])]
		return target

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
	def function_gradient_r_over_r(self,r):
		return -2.*self.gamma*self.function(r)

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
	def function_gradient(self,r):
		return -self.gamma*np.square(self.r)*np.exp(-self.gamma*self.r)
	def function_gradient_r_over_r(self,r):
		raise NotImplementedError

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
	def function_gradient_r_over_r(self,r):
		raise NotImplementedError











	
