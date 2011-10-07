from kern import kern
import numpy as np

class constant(kern):
	"""Base class for kernels which are constant in X
	to extend this calss, define a function which accepts a shape and returns an array of that shape with ones and zeros in"""
	def __init__(self,X,alpha=1.):
		kern.__init__(self)
		self.alpha = alpha
		self.set_X(X)
		self.Nparam = 1
	def set_X(self,X):
		self.shape = (X.shape[0],X.shape[0])
		self.X = X
	def gradients_X(self,target=None):
		if target is None:
			return np.zeros(self.shape+self.X.shape)

	def compute(self,target=None):
		if target is None:
			target = np.zeros(self.shape)
		if self.masked:
			target[self.mask_grid] += self.alpha*self.function((self.mask.size,self.mask.size))
		else:
			target[:,:] += self.alpha*self.function(self.shape)
		return target
	
	def gradients(self):
		return [self.function(self.shape)]
	def set_param(self,alpha):
		self.alpha = alpha
	def get_param(self):
		return self.alpha
	def get_param_names(self):
		return ['alpha']
	def diag_compute(self):
		return np.ones(self.shape[0])*self.alpha

class white(constant):
	def __init__(self,X,alpha=1.):
		constant.__init__(self,X,alpha)
	def cross_compute(self,X2,target=None):
		if target is None:
			return np.zeros((self.shape[0],X2.shape[0]))
	def function(self,shape):
		return np.eye(shape[0])

class bias(constant):
	def __init__(self,X,alpha=1.):
		constant.__init__(self,alpha)
	def function(self,shape):
		return np.ones(shape)
