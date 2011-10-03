from kern import kern

class white(kern):
	def __init__(self,X,alpha=1.):
		self.alpha = alpha
		self.set_X(X)
	def set_X(self):
		self.shape = (X.shape[0],X.shape[0])
		self.eye = np.eye(self.shape[0])
	def compute(self):
		return self.eye*self.alpha
	def gradients(self):
		return [self.eye]
	def gradients_X(self):
		pass #TODO
	def set_param(self):
		self.alpha= alpha
	def get_param(self,x):
		self.alpha = x
	def diag_compute(self):
		return np.ones(self.shape[0])*self.alpha

def bias(kern):
	def __init__(self,X,alpha=1.):
		self.alpha = alpha
		self.set_X(X)
	def set_X(self):
		self.shape = (X.shape[0],X.shape[0])
		self.ones = np.ones(self.shape)
	def compute(self):
		return self.ones*self.alpha
	def gradients(self):
		return [self.ones]
	def gradients_X(self):
		pass #TODO
	def set_param(self):
		self.alpha= alpha
	def get_param(self,x):
		self.alpha = x
	def diag_compute(self):
		return np.ones(self.shape[0])*self.alpha
	def set_param(self):
		self.alpha= alpha
	def get_param(self,x):
		self.alpha = x
