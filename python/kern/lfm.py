from kern import kern, cross_kern
from lfmUpsilon import *

class lfm(kern):
	"""The base class for the latent force model kernel"""
	def __init__(self,t,mass=1.,spring=1.,damper=1., sensitivity=1., initVal=0., rbfalpha=1., rbfgamma=1.):
		assert t.shape[1]==1,'1D inputs only for the lfm system'
		kern.__init__(self,t)
		self.mass = mass
		self.spring = spring
		self.damper = damper
		self.sensitivity = sensitivity
		self.initVal = initVal
		self.rbfalpha = rbfalpha
		self.rbfgamma = rbfgamma

	def setX(self,X):
		self.args = (X,X)

	def cross_args(self,X2):
		""" compute the arguments to the function when computing a cross variance"""
		return (self.X, X2)

	def get_param(self):
		return [self.mass, self.spring, self.damper, self.sensitivity, self.rbfgamma]

	def set_param(self, param):
		self.mass, self.spring, self.damper, self.sensitivity, self.rbfgamma = param
		self.alpha = self.damper/(2.*self.mass);
		self.omega = np.sqrt(self.spring/self.mass-self.alpha**2);
		self.gamma = self.alpha + 1j*self.omega;
		self.gamma_ = self.alpha - 1j*self.omega;
		self.zeta = self.damper/(2.*np.sqrt(self.mass*self.spring));
		self.omega_0 = np.sqrt(self.spring/self.mass);

	def get_param_names(self):
		return ['mass', 'spring', 'damper', 'sensitivity', 'rbfgamma']
	
	def function(self,t1,t2):
		"""
		TODO: There's probably a faster way to compute this, given that we're calling the same Upsilon many times in the function h. for now, this will do.
		"""
		prefactor = np.square(self.sensitivity)*np.sqrt(np.pi)*self.lengthscale/8./np.square(self.omega)
		g,g_ = self.gamma, self.gamma_
		k = h(g_,g,t1,t2) + h(g,g_,t2,t1) + h(g,g_,t1,t2) + h(g_,g,t2,t1)\
		  - h(g_,g_,t1,t2) - h(g_,g_,t2,t1) - h(g,g,t1,t2) - h(g,g,t2,t1)
		return prefactor*k

	def gradients(self,t1,t2):
		raise NotImplementedError #TODO

class lfmXlfm(cross_kern):
	def __init__(self,lfm1,lfm2):
		assert isinstance(lfm1, lfm)
		assert isinstance(lfm2, lfm)#be sure that the kernels we're crossing are both lfm
		cross_kern.__init__(self,lfm1,lfm2)
		self.mass1, self.spring1, self.damper1, self.sensitivity1, self.rbfgamma1 = lfm1.get_param()
		self.mass2, self.spring2, self.damper2, self.sensitivity2, self.rbfgamma2  = lfm2.get_param()

	def function(self,t1,t2):
		prefactor = self.sensitivity1*self.sensitivity2*np.sqrt(np.pi)*self.lengthscale/8./self.omega1/self.omega2
		g1,g1_,g2,g2_ = self.kern1.gamma, self.kern1.gamma_, self.kern2.gamma, self.kern2.gamma_
		k = h(gd_,g1,t1,t2) + h(g1,g2_,t2,t1) + h(g2,g1_,t1,t2) + h(g1_,g2,t2,t1)\
		  - h(g2_,g1_,t1,t2) - h(g1_,g2_,t2,t1) - h(g2,g1,t1,t2) - h(g1,g2,t2,t1)
		return prefactor*k

	def gradients(self,t1,t2):
		raise NotImplementedError #TODO

	def get_param(self):
		#TODO: do we need to return these here? superfluous?
		return np.array([self.mass1, self.spring1, self.damper1, self.sensitivity1, self.rbfgamma1, self.mass2, self.spring2, self.damper2, self.sensitivity2, self.rbfgamma2])

	def set_param(self,x):
		self.mass1, self.spring1, self.damper1, self.sensitivity1, self.rbfgamma1, self.mass2, self.spring2, self.damper2, self.sensitivity2, self.rbfgamma2 = x

	def get_param_names(self):
		return ['mass1', 'spring1', 'damper1', 'sensitivity1', 'rbfgamma1', 'mass2', 'spring2', 'damper2', 'sensitivity2', 'rbfgamma2']

class lfmXrbf(cross_kern):
	def __init__(self,lfmkern,rbfkern):
		cross_kern.__init__(self,lfmkern,rbfkern)

	def function(self,t1,t2):
		raise NotImplementedError #TODO

	def gradients(self,t1,t2):
		raise NotImplementedError #TODO
	
		

		
