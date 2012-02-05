from kern import kern, crosskern

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
		return [self.mass, self.spring, self.damper, self.sensitivity, self.rbfgamma

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
		raise NotImplementedError #TODO

	def gradients(self,t1,t2):
		raise NotImplementedError #TODO

class lfmXlfm(crosskern):
	def __init__(self,lfm1,lfm2):
		crosskern.__init__(self,lfm1,lfm2)

	def function(self,t1,t2):
		raise NotImplementedError #TODO

	def gradients(self,t1,t2):
		raise NotImplementedError #TODO

class lfmXrbf(crosskern):
	def __init__(self,lfmkern,rbfkern):
		crosskern.__init__(self,lfmkern,rbfkern)

	def function(self,t1,t2):
		raise NotImplementedError #TODO

	def gradients(self,t1,t2):
		raise NotImplementedError #TODO
	
		

		
