from kern import kern
import numpy as np

class constant(kern):
	"""Base class for kernels which are constant in X
	to extend this calss, define a function which accepts a shape and returns an array of that shape with ones and zeros in"""
	def __init__(self,X,alpha=1.):
		kern.__init__(self,X)
		self.alpha = np.array([alpha])
		self.expand_X(X)
	def set_X(self,X):
		self.args = ((X.shape[0],X.shape[0]),)
		self.X = X
	def extract_gradients_X(self,target=None):
		"""totally override this function, bypassing self.gradients_X. """
		return 0.
	def set_param(self,alpha):
		self.alpha = alpha
	def get_param(self):
		return self.alpha.copy()
	def get_param_names(self):
		return ['alpha']

class white(constant):
	def __init__(self,X,alpha=1.):
		constant.__init__(self,X,alpha)
		self.name='white'
	def cross_compute(self,X2,target=None):
		return np.zeros((self.shape[0],X2.shape[0]))
	def function(self,shape):
		return self.alpha*np.eye(shape[0])
	def gradients(self,shape):
		return [np.eye(shape[0])]
		

class bias(constant):
	def __init__(self,X,alpha=1.):
		constant.__init__(self,X,alpha)
		self.name='bias'
	def cross_args(self,X2):
		return ((self.shape[0],X2.shape[0]),)
	def function(self,shape):
		return self.alpha*np.ones(shape)
	def gradients(self,shape):
		return np.ones(shape)

class fixed_effect(constant):
	def __init__(self,X,K):
		constant.__init__(self,X)
		self.K = K
		self.name='fixed'
	def cross_compute(self,X2,target=None):
		return np.zeros((self.shape[0],X2.shape[0]))
	def function(self,shape):
		return self.alpha*self.K
	def gradients(self,shape):
		return [self.K]
	def compute_new(self,Xnew):
		return np.zeros((Xnew.shape[0],Xnew.shape[0]))
		

