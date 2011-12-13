import numpy as np
import pylab as pb
import sys
import kern
from simple_GP import GP 
import hgp
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
		#var = Kxx - np.dot(Kx.T,np.dot(B, Kx))
		var = Kxx - ndlutil.mdot(Kx.T,self.Kfi - ndlutil.mdot(self.Kfi,fcov,self.Kfi), Kx)
		return mu,var

	def predict_y(self,Xnew,i):
		N = self.Y.shape[1]
		B = N*self.Kyi-self.A*N**2

		fhat = np.dot(self.Li,np.dot(self.Kyi,self.Y.sum(1)))
		fcov = self.Li

		mu_f, var_f = self.predict_f(Xnew)

		Kxx = self.kerny.compute_new(Xnew)
		Kx = self.kerny.cross_compute(Xnew)

		diff = self.Y[:,i]-fhat
		
		mu = mu_f + np.dot(Kx.T,np.dot(self.Kyi,diff))
		var = var_f + Kxx - np.dot(Kx.T,np.dot(self.Kyi,Kx))
		return mu,var

	def plot(self):
		assert self.X.shape[1]==1
		ndlutil.Tango.reset()

		xmin,xmax = self.X.min(),self.X.max()
		xmin,xmax = xmin-0.2*(xmax-xmin),xmax+0.2*(xmax-xmin)
		ymin,ymax = self.Y.min(),self.Y.max()
		ymin,ymax = ymin-0.2*(ymax-ymin),ymax+0.2*(ymax-ymin)
		xx = np.linspace(xmin,xmax,100)[:,None]

		pb.subplot(1,self.Y.shape[1]+1,1)
		ndlutil.utilities.gpplot(xx,*self.predict_f(xx),alpha=0.3)
		fhat = np.dot(self.Li,np.dot(self.Kyi,self.Y.sum(1)))
		fcov = self.Li
		pb.errorbar(self.X[:,0],fhat,color='k',yerr=2*np.sqrt(np.diag(fcov)),elinewidth=2,linewidth=0)
		pb.xlim(xmin,xmax)
		pb.ylim(ymin,ymax)
		pb.xticks([])

		if isinstance(self.kernf,kern.hierarchical):
			for i,y in enumerate(self.Y.T):
				pb.subplot(1,self.Y.shape[1]+1,i+2)
				c = ndlutil.Tango.nextMedium()
				cf = ndlutil.Tango.nextLight()
				cp = ndlutil.Tango.nextDark()
				pb.plot(self.X,y[:,None],color=cp,marker='o')
				#ndlutil.utilities.gpplot(xx,*self.predict_y(xx,i),edgecol=c,fillcol=cf,alpha=0.3)
				pb.xlim(xmin,xmax)
				pb.ylim(ymin,ymax)
				pb.xticks([])
				pb.yticks([])
		else:
			for i,y in enumerate(self.Y.T):
				pb.subplot(1,self.Y.shape[1]+1,i+2)
				c = ndlutil.Tango.nextMedium()
				cf = ndlutil.Tango.nextLight()
				cp = ndlutil.Tango.nextDark()
				pb.plot(self.X,y[:,None],color=cp,marker='o')
				ndlutil.utilities.gpplot(xx,*self.predict_y(xx,i),edgecol=c,fillcol=cf,alpha=0.3)
				pb.xlim(xmin,xmax)
				pb.ylim(ymin,ymax)
				pb.xticks([])
				pb.yticks([])





		

if __name__=='__main__':

	Nd = 30
	Ng = 8
	#X = np.linspace(-3,3,Nd)[:,None]
	X = np.random.randn(Nd,1)
	X.sort(0)
	Y = np.sin(X) + 0.3*np.sin(X+5*np.random.np.random.rand(1,Ng)) + np.random.randn(Nd,Ng)*0.1
	kernf = kern.rbf(X)
	kerny = kern.rbf(X) + kern.white(X)

	m = chgp(X,Y,kernf,kerny)
	m.constrain_positive('')
	m.optimize()

	m.plot()

	#build a hgp model
	XX = np.tile(X,(Ng,1))
	YY = Y.T.flatten()[:,None]
	pd = np.vstack([np.tile([[str(i)]],(Nd,1)) for i in range(Ng)])
	m2 = hgp.hgp(XX,YY,pd)
	m2.kern.unconstrain('')
	m2.constrain_positive('')
	m2.optimize()
	pb.figure()
	m2.plot()


			
