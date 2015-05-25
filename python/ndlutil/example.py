import numpy as np
from utilities import model
import pylab as pb


class example(model):
	def __init__(self):
		model.__init__(self)
		self.set_param(np.random.randn(4))
	def get_param(self):
		return np.array([self.foo,self.bar, self.choo, self.blah])
	def set_param(self,x):
		self.foo, self.bar, self.choo, self.blah = x
	def get_param_names(self):
		return ['foo', 'bar', 'choo', 'blah']
	def log_likelihood(self):
		return -np.square(self.foo-self.bar + self.choo + self.blah)
		

