import numpy as np
import pylab as pb
import sys
import kern
from simple_GP import GP 
from ndlutil.utilities import pdinv, gpplot
import ndlutil
from scipy import optimize
import pdb


class chgp(ndlutil.model):
	def __init__(self,X,Y,kernf,kerny):
		"""
		A collapsed hierarchical GP. 

		Assuming we have several measuremnets (call them genes) sampled at the same times, 
		we can expolit this to build a GP which is more efficient for inference. 

		Arguments
		----------
		X : a Nt x D numpy array (function inputs)
		Y : a Nt x Ng numpy array, each column being a gene
		kernf: a kernel object which models connections between genes
		kerny: a kernel object which models the invididual genes

		Notes
		-----
		You can use a hgp kernel for kerny to make a double hierarchy

		"""

		assert len(X.shape)==2
		assert len(Y.shape)==2
		assert Y.shape[0] == X.shape[0]

		ndlutil.model.__init__(self)

		self.X = X
		self.Y = Y
		self.Ysum_outer = np.outer(self.Y.sum(1),self.Y.sum(1))
		self.Youter = np.dot(self.Y,self.Y.T)

		self.kernf, self.kerny = kernf, kerny
		self.NFparam = self.kernf.extract_param().size
		self.NYparam = self.kerny.extract_param().size
		self.expand_param(self.extract_param())

	def get_param(self):
		return np.hstack((self.kernf.extract_param(),self.kerny.extract_param()))

	def get_param_names(self):
		return self.kernf.extract_param_names() + self.kerny.extract_param_names()

	def set_param(self,x):
		self.kernf.expand_param(x[:self.NFparam])
		self.kerny.expand_param(x[self.NFparam:])

		self.Ky = self.kerny.compute()
		self.Kyi, self.Ky_hld = pdinv(self.Ky)
		self.Kf = self.kernf.compute() + np.eye(self.X.shape[0])*1e-6
		self.Kfi, self.Kf_hld = pdinv(self.Kf)

		N = self.Y.shape[1]
		self.Li, self.L_hld = pdinv(self.Kyi*N+self.Kfi)
		self.A = np.dot(self.Kyi,np.dot(self.Li, self.Kyi))
	
	def log_likelihood(self):
		return -0.5*self.Y.size*np.log(2.*np.pi) - self.Y.shape[1]*self.Ky_hld - self.Kf_hld - self.L_hld - 0.5*np.sum(self.Kyi*self.Youter) + 0.5*np.sum(self.A*self.Ysum_outer)

	def log_likelihood_gradients(self):
		N = self.Y.shape[1]
		B = N*self.Kyi-self.A*N**2
		dL_dKf = -0.5*self.Kfi + 0.5*np.dot(self.Kfi,np.dot(self.Li, self.Kfi)) + 0.5*np.dot(B,np.dot(self.Ysum_outer,B))/N**2
		dL_dKy = -0.5*N*self.Kyi + 0.5*np.dot(self.Kyi,np.dot(self.Li, self.Kyi))*N + 0.5*np.dot(self.Kyi,np.dot(self.Youter,self.Kyi)) +(0.5*np.dot(B,np.dot(self.Ysum_outer,B))/N**2 - 0.5*np.dot(self.Kyi,np.dot(self.Ysum_outer,self.Kyi)))/N
		return np.array(map(lambda x:np.sum(x*dL_dKf), self.kernf.extract_gradients()) + map(lambda x:np.sum(x*dL_dKy), self.kerny.extract_gradients()))

	def predict_f(self, Xnew):
		"""
		Make a prediction for the GP function. 

		"""
		N = self.Y.shape[1]
		B = N*self.Kyi-self.A*N**2

		fhat = np.dot(self.Li,np.dot(self.Kyi,self.Y.sum(1)))
		fcov = self.Li
		Kx = self.kernf.cross_compute(Xnew)
		Kxx = self.kernf.compute_new(Xnew)
		mu = np.dot(Kx.T,np.dot(self.Kfi,fhat))
		var = Kxx - np.dot(Kx.T,np.dot(B, Kx))
		return mu,var

	def predict_y(self,Xnew,i):
		N = self.Y.shape[1]
		B = N*self.Kyi-self.A*N**2

		fhat = np.dot(self.Li,np.dot(self.Kyi,self.Y.sum(1)))
		fcov = self.Li
		#TODO





		

if __name__=='__main__':

	Nd = 20
	Ng = 8
	X = np.linspace(-3,3,Nd)[:,None]
	Y = np.sin(X) + 0.3*np.sin(X+5*np.random.np.random.rand(1,Ng)) + np.random.randn(Nd,Ng)*0.01
	for y in Y.T:
		pb.plot(X,y[:,None],marker='o')
	kernf = kern.rbf(X)
	kerny = kern.rbf(X) + kern.white(X)

	m = chgp(X,Y,kernf,kerny)
	m.constrain_positive('')
	m.optimize()

	xx = np.linspace(-3.2,3.2,100)[:,None]
	ndlutil.utilities.gpplot(xx,*m.predict_f(xx))
	#pb.show()


			
