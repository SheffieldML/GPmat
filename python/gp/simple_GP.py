import numpy as np
import pylab as pb
import sys
import kern
from ndlutil import model
from ndlutil.utilities import pdinv, gpplot

class GP(model):
	def __init__(self,X,Y,kernel=None):
		model.__init__(self)
		assert len(X.shape)==2
		assert len(Y.shape)==2
		assert X.shape[0]==Y.shape[0]
		assert Y.shape[1]==1
		self.X = X
		self.Y = Y
		if kernel is None:
			self.kern = kern.rbf(self.X) + kern.white(self.X)
		else:
			self.kern = kernel
		self.Youter = np.dot(Y,Y.T)
		self.Nparam = self.kern.Nparam
	def set_param(self,p):
		self.kern.expand_param(p)
		self.K = self.kern.compute()
		self.Ki,self.hld = pdinv(self.K)
	def get_param(self):
		return self.kern.extract_param()
	def get_param_names(self):
		return self.kern.extract_param_names()
	def log_likelihood(self):
		return -self.hld -0.5*np.sum(self.Ki*self.Youter)
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
		if self.X.shape[1]==1:
			pb.figure()
			xmin,xmax = self.X.min(),self.X.max()
			xmin, xmax = xmin-0.2*(xmax-xmin), xmax+0.2*(xmax-xmin)
			Xnew = np.linspace(xmin,xmax,100)[:,None]
			m,v = self.predict(Xnew)
			gpplot(Xnew,m,v)
			pb.plot(self.X,self.Y,'kx',mew=1.5)
		elif self.X.shape[1]==2:
			xmin,xmax = self.X.min(0),self.X.max(0)
			xmin, xmax = xmin-0.2*(xmax-xmin), xmax+0.2*(xmax-xmin)
			xx,yy = np.mgrid[xmin[0]:xmax[0]:100j,xmin[1]:xmax[1]:100j]
			Xtest = np.vstack((xx.flatten(),yy.flatten())).T
			zz = self.predict(Xtest).reshape(100,100)
			pb.contour(xx,yy,zz,vmin=zz.min(),vmax=zz.max())
			pb.scatter(self.X[:,0],self.X[:,1],40,self.Y,linewidth=0)

		else:
			raise NotImplementedError, "Cannot plot GPs with more than two input dimensions"

if __name__=='__main__':
	X = np.random.randn(20,1)
	Y = np.sin(X)+np.random.randn(20,1)*0.05
	models = [GP(X,Y,k(X)) for k in kern.linear, kern.rbf, kern.mlp,kern.polynomial]#, kern.cubic_spline]
	[m.constrain_positive('') for m in models]
	[m.optimize() for m in models]
	[m.plot() for m in models]
	
	pb.show()



	
	
		
