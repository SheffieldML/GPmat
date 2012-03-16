import numpy as np
import pylab as pb
import sys
import kern
from ndlutil import model
from ndlutil.utilities import pdinv, gpplot

class vsGP(model):
	"""Variational sparse GP. Just Normal likelihoods for the moment"""
	def __init__(self,X,Y,Z, tau=0.1):
		"""
		Arguments
		---------
		X : matrix of input points (N x D)
		Y : matrix of observed points (N x 1)
		Z : matrix of inducing points (M x D)
		tau : noise precision
		foo bar
		"""

		model.__init__(self)
		self.kernX = kern.rbf(X)
		self.kernZ = kern.rbf(Z)
		self.Y = Y
		self.Youter = np.dot(Y,Y.T)

	def set_param(self,p):
		self.tau = p[0]
		self.kernX.expand_param(p[1:])
		self.kernZ.expand_param(p[1:])

		self.Knn = self.kernX.compute()
		self.Kmm = self.KernZ.compute()
		self.Kmmi, self.Kmm_hld = pdinv(self.Kmm)
		Knm = self.kernX.cross_compute(self.Z)
		#woodbury ID for inversion...
		self.Ki = np.eye(self.Y.shape[0])*self.tau - self.tau**2*np.dot(Knm,np.dot((self.Kmmi + self.tau*np.dot(Knm.T,Knm))),Knm.T)

	def get_param(self):
		return np.hstack((self.tau,self.kernX.extract_param()))
	def get_param_names(self):
		return ['noise_precision']+self.kern.extract_param_names()
	def log_likelihood(self):
		return -self.hld -0.5*np.sum(self.Ki*self.Youter) + trace_term??
	def log_likelihood_gradients(self):
		alpha = np.dot(self.Ki,self.Y.flatten())
		dL_dK = 0.5*(np.outer(alpha,alpha)-self.Ki)
		dK_dp = self.kern.extract_gradients()
		return np.array([np.sum(dK_dpi*dL_dK) for dK_dpi in dK_dp])
	def predict(self,X):
		Kx = self.kern.cross_compute(X)
		Kxx = self.kern.compute_new(X)
		mu = np.dot(np.dot(Kx.T,self.Ki),self.Y)
		var = Kxx - np.dot(np.dot(Kx.T,self.Ki),Kx)
		return mu,var
	def plot(self):
		pb.figure()
		X = self.kern.kerns[0].X
		xmin,xmax = X.min(),X.max()
		xmin, xmax = xmin-0.2*(xmax-xmin), xmax+0.2*(xmax-xmin)
		Xnew = np.linspace(xmin,xmax,100)[:,None]
		m,v = self.predict(Xnew)
		gpplot(Xnew,m,v)
		pb.plot(X,self.Y,'kx',mew=1.5)

if __name__=='__main__':
	X = np.random.randn(20,1)
	Y = np.sin(X)+np.random.randn(20,1)*0.05
	kerns = [kern.linear(X), kern.rbf(X), kern.mlp(X)]
	kerns = [k+kern.white(X) for k in kerns]
	[k.constrain_positive('') for k in kerns]
	models = [vsGP(k,Y) for k in kerns]
	[m.optimize() for m in models]
	[m.plot() for m in models]
	



	
	
		
