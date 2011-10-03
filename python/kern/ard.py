import numpy as np
from kern import kern

class ard(kern):
	def __init__(self,X,alpha=1., gamma=None):
		if gamma is None:
			gamma = np.ones(X.shape[1])
		self.alpha, self.gamma = alpha, gamma
		self.set_X(X)
	def set_X(self,X):
		self.r2 = np.square(np.sum(X[:,:,None]-X[None,:,:],-1))
		self.X = X
		self.shape = (X.shape[0],X.shape[0])
	def set_param(self,x):
		self.alpha = x[0]
		self.gamma = x[1:]
	def get_param(self):
		return np.hstack((self.alpha, self.gamma))
	def compute(self):
		return self.alpha*np.exp(-np.sum(self.r2*self.gamma[None,None,:],-1))
	def cross_compute(self,X2):
		r2 = np.square(np.sum(self.X[:,:,None]-X2[None,:,:],-1))
		return self.alpha*np.exp(-np.sum(r2*self.gamma[None,None,:],-1))
	def gradients(self):
		tmp = np.exp(-np.sum(self.r2*self.gamma[None,None,:],-1))
		return [tmp]+[-self.alpha*self.r2[:,:,i]*tmp for i in range(len(self.gamma))]
	def get_param_names(self):
		return ['alpha']+['gamma_%i'%i for i in range(len(self.gamma))]
	
	


