import numpy as np
from kern import kern

class dot(kern):
	"""
	Base class for kernels which are functions of the inner product.

	Notes
	-----
	This also allows for efficient computation: $x^{\\top}x$ is stored when set_X is called.

	Masking is handled higher up the inherritance tree, in kern

	"""
	def __init__(self,X,alpha=1.):
		kern.__init__(self,X)
		self.alpha = alpha #variance
		self.Nparam = 1
		self.set_X(X)

	def set_X(self,X):
		"""Set self.X  to the passed value and compute and store the argument(s) to the function
		Notes
		-----
		These three arguments are gramm matrices containing
		  1) $K_1[i,j] = x_i^\\top x_j$
		  2) $K_1[i,j] = x_i^\\top x_i$
		  3) $K_1[i,j] = x_j^\\top x_j$
		(actually, they aren't all full gramm matrices, but will broadcast that way!)
		  """
		self.args = (np.sum(X[:,None,:]*X[None,:,:],-1),\
				np.sum(X[:,None,:]*X[:,None,:],-1),\
				np.sum(X[None,:,:]*X[None,:,:],-1),)
		self.X = X

	def cross_args(self,X2):
		""" compute the arguments to the function when computing a cross variance"""
		return (np.sum(self.X[:,None,:]*X2[None,:,:],-1),\
				np.sum(self.X[:,None,:]*self.X[:,None,:],-1),\
				np.sum(X2[None,:,:]*X2[None,:,:],-1),)

	def function(self):
		raise NotImplementedError
	def gradients(self):
		raise NotImplementedError
	def gradients_X(self):
		raise NotImplementedError, "TODO"
	def get_param(self):
		return self.alpha
	def set_param(self,x):
		self.alpha = x
	def get_param_names(self):
		return ['alpha']

class linear(dot):
	"""
	A linear Kernel
	Notes
	-----
	$$
	k(x,x') = \\alpha x^\\top x'
	$$
	"""
	def __init__(self,X,alpha=1.):
		dot.__init__(self,X,alpha)
	def function(self,xTx,x1x1,x2x2):
		return self.alpha*xTx
	def gradients(self,xTx,x1x1,x2x2):
		return [xTx]

class polynomial(dot):
	"""A polynomial kernel with a fixed order"""
	def __init__(self,X,alpha=1.,order=3.):
		dot.__init__(self,X,alpha)
		self.order=order
	def function(self,xTx,x1x1,x2x2):
		return self.alpha*(1.+xTx)**self.order
	def gradients(self,xTx,x1x1,x2x2):
		return [(1.+xTx)**self.order]

class mlp(dot):
	"""
	A mlp kernel. This has an extra parameter over the linear kernel.
	
	Notes
	-----
	It's possible to have an ARD implementation of this kernel,
	but I'm going to implement ARD as a separate function, applicable to all any mutlivariate kernel. TODO!
	"""
	def __init__(self,X,alpha=1.,gamma=3.):
		dot.__init__(self,X,alpha)
		self.gamma=gamma
		self.Nparam = 2
	def function(self,xTx,x1x1,x2x2):
		return self.alpha*2./np.pi*np.arcsin(self.gamma*xTx/np.sqrt((1.+self.gamma*x1x1)*(1.+self.gamma*x2x2)))
	def gradients(self,xTx,x1x1,x2x2):
		"""Credit to wolfram alpha :)"""
		denom2 = (1.+self.gamma*x1x1)*(1.+self.gamma*x2x2)
		return [2./np.pi*np.arcsin(self.gamma*xTx/np.sqrt(denom2)),\
			self.alpha/np.pi*xTx*(self.gamma*(x1x1+x2x2)+2.)/np.power(denom2,3./2.)/np.sqrt(1.-np.square(self.gamma*xTx)/denom2)]
	def get_param(self):
		return np.array([self.alpha,self.gamma]).flatten()
	def set_param(self,x):
		self.alpha,self.gamma = x
	def get_param_names(self):
		return ['alpha', 'gamma']



	

