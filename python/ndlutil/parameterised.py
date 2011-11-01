import numpy as np
from scipy import linalg, optimize
import pylab as pb
import Tango
import sys
import re
import numdifftools as ndt
import pdb
import cPickle
from utilities import sigmoid


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
	
	def pickle(self,filename,protocol=-1):
		f = file(filename,'w')
		cPickle.dump(self,f,protocol)
		f.close()

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
		"""Return a np array of all the constrained indices"""
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
			if ((xx<=lower)|(xx>=upper)) & (i in matches):
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
		[np.put(x,i,np.log(np.clip(x[i]-l,1e-10,np.inf)/np.clip(h-x[i],1e-10,np.inf))) for i,l,h in zip(self.constrained_bounded_indices, self.constrained_bounded_lowers, self.constrained_bounded_uppers)]

		to_remove = self.constrained_fixed_indices+[t[1:] for t in self.tied_indices]
		if len(to_remove):
			return np.delete(x,np.hstack(to_remove))
		else:
			return x


	def expand_param(self,x):
		""" takes the vector x, which is then modified (by untying, reparameterising or inserting fixed values), and then call self.set_param"""
		#work out how many places are fixed, and where they are. tricky logic!
		Nfix_places = 0.
		if len(self.tied_indices):
			Nfix_places += np.hstack(self.tied_indices).size-len(self.tied_indices)
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



