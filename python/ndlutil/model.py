import numpy as np
from scipy import linalg, optimize
import pylab as pb
import sys
import re
import numdifftools as ndt
import pdb
from parameterised import parameterised
import prior
from utilities import sigmoid

class model(parameterised):
	def __init__(self):
		parameterised.__init__(self)
	def get_param(self):
		raise NotImplementedError, "this needs to be implemented to utilise the model class"
	def set_param(self,x):
		raise NotImplementedError, "this needs to be implemented to utilise the model class"
	def log_likelihood(self):
		raise NotImplementedError, "this needs to be implemented to utilise the model class"
	def log_likelihood_gradients(self):
		raise NotImplementedError, "this needs to be implemented to utilise the model class"

	def set_prior(self,which,what):
		"""sets priors on the model parameters.

		Arguments
		---------
		which -- string, regexp, or integer array
		what -- instance of a prior class
		"""
		if not hasattr(self,'priors'):
			self.priors = [None for i in range(self.get_param().size)]

		which = self.grep_param_names(which)

		#priors on positive variables
		if isinstance(what, (prior.gamma, prior.log_Gaussian)):
			assert not np.any(which[:,None]==self.constrained_negative_indices), "constraint and prior incompatible"
			unconst = np.setdiff1d(which, self.constrained_positive_indices)
			if len(unconst):
				print "Warning: constraining parameters to be positive:"
				print '\n'.join([n for i,n in enumerate(self.get_param_names()) if i in unconst])
				print '\n'
				self.constrain_positive(unconst)


		#store the prior in a local list
		for w in which:
			self.priors[w] = what


	def log_prior(self):
		"""evaluate the prior"""
		if not hasattr(self, 'priors'):
			return 0.
		return np.sum([p.lnpdf(x) for p, x in zip(self.priors,self.get_param()) if p is not None])
	
	def log_prior_gradients(self):
		"""evaluate the gradients of the priors"""
		x = self.get_param()
		ret = np.zeros(x.size)
		if hasattr(self, 'priors'):
			[np.put(ret,i,p.lnpdf_grad(xx)) for i,(p,xx) in enumerate(zip(self.priors,x)) if not p is None]
		return ret

	def extract_gradients(self):
		"""use self.log_likelihood_gradients and self.prior_gradients to get the gradients of the model.
		Adjust the gradient for constraints and ties, return."""
		g = self.log_likelihood_gradients() + self.log_prior_gradients()
		x = self.get_param()
		g[self.constrained_positive_indices] = g[self.constrained_positive_indices]*x[self.constrained_positive_indices]
		g[self.constrained_negative_indices] = g[self.constrained_negative_indices]*x[self.constrained_negative_indices]
		[np.put(g,i,g[i]*(x[i]-l)*(h-x[i])/(h-l)) for i,l,h in zip(self.constrained_bounded_indices, self.constrained_bounded_lowers, self.constrained_bounded_uppers)]
		[np.put(g,i,v) for i,v in [(t[0],np.sum(g[t])) for t in self.tied_indices]]
		if len(self.tied_indices):
			to_remove = np.hstack((self.constrained_fixed_indices,np.hstack([t[1:] for t in self.tied_indices])))
		else:
			to_remove=self.constrained_fixed_indices
			
		return np.delete(g,to_remove)

	def randomize(self):
		"""
		Randomize the model. 
		Make this draw from the prior if one exists, else draw from N(0,1)
		"""
		#first take care of all parameters (from N(0,1))
		x = self.extract_param()
		x = np.random.randn(x.size)
		self.expand_param(x)
		#now draw from prior where possible
		if hasattr(self,'priors'):
			x = self.get_param()
			[np.put(x,i,p.rvs(1)) for i,p in enumerate(self.priors) if not p is None]
			self.set_param(x)


	def optimize_restarts(self, Nrestarts=10, compare_Laplace=False, **kwargs):
		"""
		Perform random restarts of the model. Include the current point in the comparison.
		passes **kwargs to the optimizer
		"""
		scores = [self.log_likelihood()]
		params = [self.extract_param()]
		for i in range(Nrestarts):
			self.randomize()
			try:
				self.optimize(**kwargs)
			except np.linalg.LinAlgErr:
				scores.append(-np.inf)
				params.append('')
				continue
			if compare_Laplace:
				self.optimize(ftol=1e-9)#need more numerical stability for good laplace approximation
				scores.append(self.Laplace_evidence())
			else:
				scores.append(self.log_likelihood())
			params.append(self.extract_param())
		i = np.argmax(scores)
		self.expand_param(params[i])

	def optimize(self,**kwargs):
		"""
		Optimize the model using self.log_likelihood and self.log_likelihood_gradient, as well as self.priors.
		Starts from the current condition and runs scipy.optimize.fmin_tnc with the passed kwargs. 
		"""
		def f_fp(x):
			try:
				self.expand_param(x)
				return -self.log_likelihood(),-self.extract_gradients()
			except:
				return np.nan, x*np.nan
		start = self.extract_param()
		opt = optimize.fmin_tnc(f_fp,start,**kwargs)[0]
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

	def checkgrad(self, include_priors=False, step=1e-6, tolerance = 1e-3, *args):
		"""check the gradient of the model by comparing to a numerical estimate. 
		If the overall gradient fails, invividual components are tested numerically#
		
		TODO: check the gradients of priors"""
	
		x = self.extract_param().copy()

		#choose a random direction to step in:
		dx = step*np.sign(np.random.uniform(-1,1,x.size))
	
		#evaulate around the point x
		self.expand_param(x+dx)
		f1,g1 = self.log_likelihood() + self.log_prior(), self.extract_gradients()
		self.expand_param(x-dx)
		f2,g2 = self.log_likelihood() + self.log_prior(), self.extract_gradients()
		self.expand_param(x)
		gradient = self.extract_gradients()
		
		numerical_gradient = (f1-f2)/(2*dx)
		ratio = (f1-f2)/(2*np.dot(dx,gradient))
		print "gradient = ",gradient
		print "numerical gradient = ",numerical_gradient
		print "ratio = ", ratio, '\n'
		sys.stdout.flush()
		
		if np.abs(1.-ratio)>tolerance:
			print "Ratio far from unity. Testing individual gradients\n"
			try:
				names = self.extract_param_names()
			except NotImplementedError:
				names = ['Variable %i'%i for i in range(len(x))]
			for i in range(len(x)):
				xx = x.copy()
				xx[i] += step
				self.expand_param(xx)
				f1,g1 = self.log_likelihood() + self.log_prior(), self.extract_gradients()[i]
				xx[i] -= 2.*step
				self.expand_param(xx)
				f2,g2 = self.log_likelihood() + self.log_prior(), self.extract_gradients()[i]
				self.expand_param(x)
				gradient = self.extract_gradients()[i]

			
				numerical_gradient = (f1-f2)/(2*step)
				print names[i]
				ratio = (f1-f2)/(2*step*gradient)
				difference = np.abs((f1-f2)/2/step - gradient)
				print "ratio = ",ratio
				print "difference = ",difference,'\n'
				sys.stdout.flush()
			return False
		return True
	
