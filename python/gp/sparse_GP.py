import numpy as np
import pylab as pb
import sys
import kern
from ndlutil import model
from ndlutil.utilities import pdinv, gpplot

class vsGP(model):
	def __init__(self,X,Y,M=10, tau=50.):
		model.__init__(self)
		assert len(X.shape)==2
		assert len(Y.shape)==2
		assert X.shape[0]==Y.shape[0]
		assert Y.shape[1]==1 #TODO: multiple independent ops
		self.X = X
		self.N, self.D = self.X.shape
		self.M = M
		self.Z = self.X[np.random.permutation(self.N)[:self.M]]
		self.Y = Y
		self.tau = tau
		self.kernN = kern.rbf(self.X)
		self.kernM = kern.rbf(self.Z)
		self.Youter = np.dot(Y,Y.T)

	def set_param(self,p):
		self.tau = p[0]
		self.kernN.expand_param(p[1:])
		self.kernM.expand_param(p[1:])
		self.Knn_diag = np.diag(self.kernN.compute()) # TODO: make the kern toolbox compute diags
		self.Kmm = self.kernM.compute()
		self.Kmmi, self.Kmm_hld = pdinv(self.Kmm)
		self.Knm = self.kernN.cross_compute(self.Z)
		self.KmnKnm = np.dot(self.Knm.T, self.Knm)
		self.KnmKmmiKmn = np.dot(self.Knm,np.dot(self.Kmmi, self.Knm.T))
		self.Woodbury_inv, self.Woodbury_hld = pdinv(self.tau*self.KmnKnm + self.Kmm)
		self.Qi = np.dot(self.Knm, np.dot(self.Woodbury_inv,self.Knm.T))
		self.hld = -0.5*self.N*np.log(self.tau) +self.Kmm_hld + self.Woodbury_hld

	def get_param(self):
		return np.hstack([self.tau,self.kernM.extract_param()])
	def get_param_names(self):
		return ['tau']+self.kernM.extract_param_names()
	def log_likelihood(self):
		return -0.5*self.Y.size*np.log(2.*np.pi) -self.hld -0.5*np.sum(self.Qi*self.Youter) -0.5*self.tau*(self.Knn_diag.sum() + np.trace(self.KnmKmmiKmn))
	def log_likelihood_gradients(self):
		pass #return np.array([np.sum(dK_dpi*dL_dK) for dK_dpi in dK_dp])
	def predict(self,X):
		u_mu = self.tau*np.dot(self.Kmm,np.dot(self.Woodbury_inv,np.dot(self.Knm.T,self.Y))) # mean of the latent variables f_m (or u)
		Kx = self.kernM.cross_compute(X)
		Kxx = self.kernM.compute_new(X)
		mu = np.dot(np.dot(Kx.T,self.Kmmi),u_mu)
		var = Kxx - np.dot(np.dot(Kx.T,self.Kmmi-self.Woodbury_inv),Kx) 
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
			pb.plot(self.Z,self.Z*0+pb.ylim()[0],'k|',mew=1.5,markersize=12)
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

