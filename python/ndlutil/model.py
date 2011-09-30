import numpy as np
from scipy import linalg, optimize
import pylab as pb
import Tango
import sys
import re
import numdifftools as ndt
import pdb
import cPickle
from parameterised import parameterised


class model(parameterised):
	def __init__(self):
		parameterised.__init__(self)

	def set_prior(self,which,what):
		"""sets priors on the model parameters.

		Arguments
		---------
		which -- string, regexp, or integer array
		what -- instance of prior class
		"""
		pass

	def log_prior(self):
		"""evaluate the prior"""
		raise NotImplementedError, "TODO"
		return self.prior.lnpdf(self.get_param())
	
	def log_prior_gradients(self):
		"""evaluate the gradients of the prior"""
		raise NotImplementedError, "TODO"
		return self.prior.lnpdf_grad(self.get_param())

	def randomize(self):
		x = self.extract_param()
		self.expand_param(np.random.randn(x.size))

	def maximum_likelihood_restarts(self, Nrestarts=10, compare_Laplace=False, **kwargs):
		D = self.extract_param().size
		scores = []
		params = []
		for i in range(Nrestarts):
			self.expand_param(np.random.randn(D))
			self.maximum_likelihood(**kwargs)
			if compare_Laplace:
				self.maximum_likelihood(ftol=1e-9)#need more numerical stability for good laplace approximation
				scores.append(self.Laplace_evidence())
			else:
				scores.append(self.log_likelihood())
			params.append(self.extract_param())
		i = np.argmax(scores)
		self.expand_param(params[i])

	def maximum_likelihood(self,**kwargs):
		def f_fp(x):
			self.expand_param(x)
			return -self.log_likelihood(),-self.extract_gradients()
		start = self.extract_param()
		opt = optimize.fmin_tnc(f_fp,start,**kwargs)[0]
		self.expand_param(opt)

	def maximum_aposteriori(self):
		def f_fp(x):
			self.expand_param(x)
			return -self.log_likelihood() -self.log_prior() ,-self.log_likelihood_gradients() -self.log_prior_gradients()
		start = self.extract_param()
		opt = optimize.fmin_tnc(f_fp,start)[0]
		self.expand_param(opt)


	def Laplace_covariance(self):
		"""return the covariance matric of a Laplace approximatino at the current (stationary) point"""
		#TODO add in the prior contributions for MAP estimation
		#TODO fix the hessian for tied, constrained and fixed components
		if hasattr(self, 'log_likelihood_hessian'):
			A = -self.log_likelihood_hessian()

		else:
			print "numerically calculating hessian. please be patient!"
			x = self.get_param()
			def f(x):
				self.set_param(x)
				return self.log_likelihood()
			h = ndt.Hessian(f)
			A = -h(x)
			self.set_param(x)
		# check for almost zero components on the diagonal which screw up the cholesky
		aa = np.nonzero((np.diag(A)<1e-6) & (np.diag(A)>0.))[0]
		A[aa,aa] = 0.
		return A

	def Laplace_evidence(self):
		"""Returns an estiamte of the model evidence based on the Laplace approximation. 
		Uses a numerical estimate of the hessian if none is available analytically"""
		A = self.Laplace_covariance()
		try:
			hld = np.sum(np.log(np.diag(jitchol(A)[0])))
		except:
			return np.nan
		return 0.5*self.get_param().size*np.log(2*np.pi) + self.log_likelihood() - hld

	def checkgrad(self,step=1e-6, tolerance = 1e-3, *args):
		"""check the gradient of the model by comparing to a numerical estimate. 
		If the overall gradient fails, invividual components are tested numerically"""
	
		x = self.get_param().copy()

		#choose a random direction to step in:
		dx = step*np.sign(np.random.uniform(-1,1,x.size))
	
		#evaulate around the point x
		self.set_param(x+dx)
		f1,g1 = self.log_likelihood(), self.log_likelihood_gradients()
		self.set_param(x-dx)
		f2,g2 = self.log_likelihood(), self.log_likelihood_gradients()
		self.set_param(x)
		gradient = self.log_likelihood_gradients()
		
		numerical_gradient = (f1-f2)/(2*dx)
		ratio = (f1-f2)/(2*np.dot(dx,gradient))
		print "gradient = ",gradient
		print "numerical gradient = ",numerical_gradient
		print "ratio = ", ratio, '\n'
		sys.stdout.flush()
		
		if np.abs(1.-ratio)>tolerance:
			print "Ratio far from unity. Testing individual gradients"
			for i in range(len(x)):
				dx = np.zeros(x.shape)
				dx[i] = step*np.sign(np.random.uniform(-1,1,x[i].shape))
				
				self.set_param(x+dx)
				f1,g1 = self.log_likelihood(), self.log_likelihood_gradients()
				self.set_param(x-dx)
				f2,g2 = self.log_likelihood(), self.log_likelihood_gradients()
				self.set_param(x)
				gradient = self.log_likelihood_gradients()

			
				numerical_gradient = (f1-f2)/(2*dx)
				print i,"th element"
				#print "gradient = ",gradient
				#print "numerical gradient = ",numerical_gradient
				ratio = (f1-f2)/(2*np.dot(dx,gradient))
				print "ratio = ",ratio,'\n'
				sys.stdout.flush()
	
