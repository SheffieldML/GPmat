import numpy as np
from scipy import linalg, optimize
import pylab as pb
import Tango
import sys
import re
import numdifftools as ndt
import pdb

class parameterised:
	def __init__(self):
		"""
		This is the base class for model and kernel. Mostly just handles tieing and constraining of parameters
		"""
		self.tied_indices = []
		self.constrained_fixed_indices = []
		self.constrained_fixed_values = []
		self.constrained_positive_indices = np.empty(shape=(0,),dtype=np.int64)
		self.constrained_negative_indices = np.empty(shape=(0,),dtype=np.int64)
		self.constrained_bounded_indices = []
		self.constrained_bounded_uppers = []
		self.constrained_bounded_lowers = []

	def tie_param(self, which):
		matches = self.grep_param_names(which)
		assert matches.size > 1, "need at least two things to tie together"
		if len(self.tied_indices):
			assert not np.any(matches[:,None]==np.hstack(self.tied_indices)), "Some indices are already tied!"
		self.tied_indices.append(matches)
		#TODO only one of the priors will be evaluated. Give a warning message if the priors are not identical
		if hasattr(self,'prior'):
			pass

		self.expand_param(self.extract_param())# sets tied parameters to single value
		

	def all_constrained_indices(self):
		if len(self.constrained_bounded_indices):
			return np.hstack([self.constrained_positive_indices, self.constrained_negative_indices, np.hstack(self.constrained_bounded_indices)])
		else:
			return np.hstack([self.constrained_positive_indices, self.constrained_negative_indices])

	def grep_param_names(self, expr):
		"""
		Arguments
		---------
		expr -- can be a regular expression object or a string to be turned into regular expression object. 

		Returns
		-------
		the indices of self.get_param_names which match the regular expression. 

		Notes
		-----
		Other objects are passed through - i.e. integers which were'nt meant for grepping
		"""

		if type(expr) is str:
			expr = re.compile(expr)
			return np.nonzero([expr.search(name) for name in self.get_param_names()])[0]
		elif type(expr) is re._pattern_type:
			return np.nonzero([expr.search(name) for name in self.get_param_names()])[0]
		else:
			return expr


	def constrain_positive(self, which):
		"""
		Set positive constraints. 

		Arguments
		---------
		which -- np.array(dtype=int), or regular expression object or string
		"""
		matches = self.grep_param_names(which)
		assert not np.any(matches[:,None]==self.all_constrained_indices()), "Some indices are already constrained"
		self.constrained_positive_indices = np.hstack((self.constrained_positive_indices, matches))
		#check to ensure constraint is in place
		x = self.get_param()
		for i,xx in enumerate(x):
			if (xx<0) & (i in matches):
				x[i] = np.exp(xx)
		self.set_param(x)


	def constrain_negative(self,which):
		"""
		Set negative constraints. 

		Arguments
		---------
		which -- np.array(dtype=int), or regular expression object or string
		"""
		matches = self.grep_param_names(which)
		assert not np.any(matches[:,None]==self.all_constrained_indices()), "Some indices are already constrained"
		self.constrained_negative_indices = np.hstack(self.constrained_positive_indices, matches)
		#check to ensure constraint is in place
		x = self.get_param()
		for i,xx in enumerate(x):
			if (xx>0) & (i in matches):
				x[i] = -np.exp(xx)
		self.set_param(x)



	def constrain_bounded(self, which, lower, upper):
		"""Set bounded constraints. 

		Arguments
		---------
		which -- np.array(dtype=int), or regular expression object or string
		upper -- (float) the upper bound on the constraint
		lower -- (float) the lower bound on the constraint
		"""
		matches = self.grep_param_names(which)
		assert not np.any(matches[:,None]==self.all_constrained_indices()), "Some indices are already constrained"
		assert lower < upper, "lower bound must be smaller than upper bound!"
		self.constrained_bounded_indices.append(matches)
		self.constrained_bounded_uppers.append(upper)
		self.constrained_bounded_lowers.append(lower)
		#check to ensure constraint is in place
		x = self.get_param()
		for i,xx in enumerate(x):
			if ((xx<lower)|(xx>upper)) & (i in matches):
				x[i] = sigmoid(xx)*(upper-lower) + lower
		self.set_param(x)


	def constrain_fixed(self, which, value):
		"""
		Arguments
		---------
		which -- np.array(dtype=int), or regular expression object or string
		value -- a float to fix the matched values to

		Notes
		-----
		Fixing a parameter which is tied to another, or constrained in some way will result in an error. 
		To fix multiple parameters to the same value, simply pass a regular expression which matches both parameter names, or pass both of the indexes
		"""
		matches = self.grep_param_names(which)
		assert not np.any(matches[:,None]==self.all_constrained_indices()), "Some indices are already constrained"
		self.constrained_fixed_indices.append(matches)
		self.constrained_fixed_values.append(value)


	def extract_param(self):
		"""use self.get_param to get the 'true' parameters of the model, which are then tied, constrained and fixed"""
		x = self.get_param()
		x[self.constrained_positive_indices] = np.log(x[self.constrained_positive_indices])
		x[self.constrained_negative_indices] = np.log(-x[self.constrained_negative_indices])
		[np.put(x,i,np.log((x[i]-l)/(h-x[i]))) for i,l,h in zip(self.constrained_bounded_indices, self.constrained_bounded_lowers, self.constrained_bounded_uppers)]
		if len(self.tied_indices):
			to_remove = np.hstack((self.constrained_fixed_indices,np.hstack([t[1:] for t in self.tied_indices])))
		else:
			to_remove=self.constrained_fixed_indices
		return np.delete(x,to_remove)

	def extract_gradients(self):
		"""use self.log_likelihood_gradients and self.prior_gradients to get the gradients of the model.
		Adjust the gradient for constraints and ties, return."""
		#TODO: prior gradients
		g = self.log_likelihood_gradients()
		x = self.get_param()
		g[self.constrained_positive_indices] = g[self.constrained_positive_indices]*x[self.constrained_positive_indices]
		g[self.constrained_negative_indices] = g[self.constrained_negative_indices]*x[self.constrained_negative_indices]
		[np.put(g,i,g[i]*(1.-sigmoid(xx[i]))*sigmoid(xx[i])*(high-low)) for i,low,high in zip(self.constrained_bounded_indices, self.constrained_bounded_lowers, self.constrained_bounded_uppers)]
		[np.put(g,i,v) for i,v in [(t[0],np.sum(g[t[1:]])) for t in self.tied_indices]]
		if len(self.tied_indices):
			to_remove = np.hstack((self.constrained_fixed_indices,np.hstack([t[1:] for t in self.tied_indices])))
		else:
			to_remove=self.constrained_fixed_indices
			
		return np.delete(g,to_remove)



	def expand_param(self,x):
		""" takes the vector x, which is then modified (by untying, reparameterising or inserting fixed values), and then call self.set_param"""
		#work out how many places are fixed, and where they are. tricky logic!
		Nfix_places = 0.
		if len(self.tied_indices):
			Nfix_places += np.hstack(self.tied_indices)-len(self.tied_indices)
		if len(self.constrained_fixed_indices):
			Nfix_places += np.hstack(self.constrained_fixed_indices).size
		if Nfix_places:
			fix_places = np.hstack(self.constrained_fixed_indices+[t[1:] for t in self.tied_indices])
		else:
			fix_places = []
		free_places = np.setdiff1d(np.arange(Nfix_places+x.size,dtype=np.int),fix_places)
		#put the models values in the vector xx
		xx = np.zeros(Nfix_places+free_places.size,dtype=np.float64)
		xx[free_places] = x
		[np.put(xx,i,v) for i,v in zip(self.constrained_fixed_indices, self.constrained_fixed_values)]
		[np.put(xx,i,v) for i,v in [(t[1:],xx[t[0]]) for t in self.tied_indices] ]
		xx[self.constrained_positive_indices] = np.exp(xx[self.constrained_positive_indices])
		xx[self.constrained_negative_indices] = -np.exp(xx[self.constrained_negative_indices])
		[np.put(xx,i,low+sigmoid(xx[i])*(high-low)) for i,low,high in zip(self.constrained_bounded_indices, self.constrained_bounded_lowers, self.constrained_bounded_uppers)]
		self.set_param(xx)

	def extract_param_names(self):
		"""Resulting parameter names after tieing."""
		n =  self.get_param_names()
		if len(self.tied_indices):
			for t in self.tied_indices:
				n[t[0]] = "<tie>".join([n[tt] for tt in t])
			remove = np.hstack([t[1:] for t in self.tied_indices])
			n = [nn for i,nn in enumerate(n) if not i in remove]
		return n

	def __str__(self):
		"""Return a string describing the parameter names and their ties and constraints"""
		#TODO



class model(parameterised):
	def __init__(self):
		parameterised.__init__(self)

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
		hld = np.sum(np.log(np.diag(jitchol(A)[0])))
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
		

class Metropolis_Hastings:
	def __init__(self,model,cov=None):
		"""Metropolis Hastings, with tunings according to Gelman et al. """
		self.model = model
		current = self.model.extract_param()
		self.D = current.size
		self.chains = []
		if cov is None:
			self.cov = m.Laplace_covariance()
		else:
			self.cov = cov
		self.scale = 2.4/np.sqrt(self.D)
		self.new_chain(current)

	def new_chain(self, start=None):
		self.chains.append([])
		if start is None:
			self.model.randomise()
		else:
			self.model.expand_param(start)



	def sample(self, Ntotal, Nburn, Nthin, tune=True, tune_throughout=False, tune_interval=100):
		current = self.model.extract_param()
		fcurrent = self.model.log_likelihood() # TODO priors
		accepted = np.zeros(Ntotal,dtype=np.bool)
		for it in range(Ntotal):
			print "sample %d of %d\r"%(it,Ntotal)
			sys.stdout.flush()
			prop = current + np.random.multivariate_normal(np.zeros(self.D),self.cov*self.scale)
			self.model.expand_param(prop)
			fprop = self.model.log_likelihood()

			if fprop>fcurrent:#sample accepted, going 'uphill'
				accepted[it] = True
				current = prop
				fcurrent = fprop
			else:
				u = np.random.rand()
				if np.exp(fprop-fcurrent)>u:#sample accepted downhill
					accepted[it] = True
					current = prop
					fcurrent = fprop

			#store current value
			if (it > Nburn) & ((it%Nthin)==0):
				self.chains[-1].append(current)

			#tuning!
			if ((it%tune_interval)==0) & tune & ((it<Nburn) | (tune_throughout)):
				pc = np.mean(accepted[it-tune_interval:it])
				if pc > .25:
					self.scale *= 1.1
				if pc < .15:
					self.scale /= 1.1

	def predict(self,function,args):
		"""Make a prediction for the function, to which we will pass the additional arguments"""
		param = self.model.get_param()
		fs = []
		for p in self.chain:
			self.model.set_param(p)
			fs.append(function(*args))
		self.model.set_param(param)# reset model to starting state
		return fs






def gpplot(x,mu,var,edgecol='k',fillcol=Tango.coloursHex['Aluminium3'],**kwargs):
	mu = mu.flatten()
	x = x.flatten()
	pb.plot(x,mu,color=edgecol,linewidth=2)
	if len(var.shape)>1:
		err = 2*np.sqrt(np.diag(var))
	else:
		err = 2*np.sqrt(var)

	kwargs['linewidth']=0
	pb.fill(np.hstack((x,x[::-1])),np.hstack((mu+err,mu[::-1]-err[::-1])),color=fillcol,**kwargs)
	#this is the edge:
	pb.plot(x,mu+err,color=edgecol,linewidth=0.2)
	pb.plot(x,mu-err,color=edgecol,linewidth=0.2)



def jitchol(A,maxtries=5):
	"""
	Arguments
	---------
	A : An almost pd square matrix

	Returns
	-------
	cho_factor(K)

	Notes
	-----
	Adds jitter to K, to enforce positive-definiteness
	"""
	
	try:
		return linalg.cho_factor(A)
	except:
		diagA = np.diag(A)
		if np.any(diagA<0.):
			raise linalg.LinAlgError, "not pd: negative diagonal elements"
		jitter= diagA.mean()*1e-6
		for i in range(1,maxtries+1):
			try:
				return linalg.cho_factor(A+np.eye(A.shape[0])*jitter)
			except:
				jitter *= 10
		raise linalg.LinAlgError,"not positive definite, even with jitter."


def pdinv(A):
	"""
	Arguments
	---------
	A : A DxD pd numpy array

	Returns
	-------
	inv : the inverse of A
	hld: 0.5* the log of the determinant of A
	"""
	L = jitchol(A)
	hld = np.sum(np.log(np.diag(L[0]))) 
	inv = linalg.flapack.dpotri(L[0],True)[0] 
	inv = np.triu(inv)+np.triu(inv,1).T 
	return inv, hld


def multiple_pdinv(A):
	"""
	Arguments
	---------
	A : A DxDxN numpy array (each A[:,:,i] is pd)

	Returns
	-------
	invs : the inverses of A
	hld: 0.5* the log of the determinants of A
	"""
	N = A.shape[-1]
	chols = [jitchol(A[:,:,i]) for i in range(N)]
	halflogdets = [np.sum(np.log(np.diag(L[0]))) for L in chols]
	invs = [linalg.flapack.dpotri(L[0],True)[0] for L in chols]
	invs = [np.triu(I)+np.triu(I,1).T for I in invs]
	return np.dstack(invs),np.array(halflogdets)

def sigmoid(x):
	return 1./(1.+np.exp(-x))

def softmax(x):
	ex = np.exp(x)
	return ex/ex.sum(1)[:,np.newaxis]

def single_softmax(x):
	ex = np.exp(x)
	return ex/ex.sum()


