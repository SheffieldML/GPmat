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

	Masking is handled higher up the inherritance tree, in kern

	"""
	def __init__(self,X,alpha=1., gamma=1.):
		kern.__init__(self,X)
		self.alpha = alpha #variance
		self.gamma = gamma #parameter of function, to be defined by inherreted class ( usually a lengthscale)
		self.Nparam = 2
		self.set_X(X)

	def set_X(self,X):
		"""Set self.X  to the passed value and compute and store the argument(s) to the function """
		self.args = (np.sqrt(np.sum(np.square(X[:,None,:]-X[None,:,:]),-1)),)
		self.X = X

	def cross_args(self,X2):
		""" compute the arguments to the function when computing a cross variance"""
		return (np.sqrt(np.sum(np.square(self.X[:,None,:]-X2[None,:,:]),-1)),)

	def function(self):
		raise NotImplementedError
	def gradients(self):
		raise NotImplementedError
	def gradients_r_over_r(self):
		raise NotImplementedError

	def gradients_X(self):
		"""
		Relies on the implementation of gradients_r_over_r:
		$$
		\\frac{\\partial f(r_{i,j}}{\\partial r_{i,j} \\times \\frac{1}{r_{i,j}}
		$$
		which arises from applying the chain rule to the derivative. Defining it this way helps prevent divide-by-zero errors

		"""
		return self.gradients_r_over_r(*self.args)[:,:,None]*(self.X[None,:,:]-self.X[:,None,:])

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
		return self.alpha*np.exp(-self.gamma*np.square(r))
	def gradients(self,r):
		f = np.exp(-self.gamma*np.square(r))
		return [f,-self.alpha*np.square(r)*f]
	def gradients_r_over_r(self,r):
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
		return (1. + self.gamma*r)*np.exp(-self.gamma*r)
	def gradients(self,r):
		f = self.function(r)
		return [f,-self.gamma*np.square(r)*f]
	def gradients_r_over_r(self,r):
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
		return self.alpha*np.exp(-self.gamma*r)
	def gradients(self):
		f = np.exp(-self.gamma*r)
		return [f, -self.gamma*r*f]
	def gradients_r_over_r(self,r):
		raise NotImplementedError











	
