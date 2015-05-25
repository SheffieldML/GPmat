import numpy as np
from scipy import linalg, optimize
import pylab as pb
import Tango
import sys
import re
import numdifftools as ndt
import pdb
import cPickle



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
				print 'Warning: adding jitter of '+str(jitter)
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


