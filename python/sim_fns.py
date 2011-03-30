import numpy as np
import pylab as pb
import pdb
from scipy.special import erf, erfc
from erfcx import erfcx


def lnDiffErfs(x1,x2):
	#function [v, signs] = lnDiffErfs(x1, x2),

	# LNDIFFERFS Helper function for computing the log of difference
	#   of two erfs.
	# FORMAT
	# DESC computes the log of the difference of two erfs in a numerically stable manner.
	# ARG x1 : argument of the positive erf
	# ARG x2 : argument of the negative erf
	# RETURN v : log(abs(erf(x1) - erf(x2)))
	# RETURN s : sign(erf(x1) - erf(x2))
	#
	# FORMAT
	# DESC computes the log of the difference of two erfs in a numerically stable manner.
	# ARG x1 : argument of the positive erf
	# ARG x2 : argument of the negative erf
	# RETURN v : log(erf(x1) - erf(x2))     (Can be complex)
	#
	# COPYRIGHT : Antti Honkela, 2007, 2008
	#
	# MODIFICATIONS : David Luengo, 2009
	#
	# MODIFICATIONS : James Hensman, 2011
	#
	# SEEALSO : gradLnDiffErfs


	x1 = np.atleast_1d(x1)
	x2 = np.atleast_1d(x2)

	x1 = x1.real
	x2 = x2.real

	#TODO make this work in the numpy style, with broadcasting and everything
	assert((x1.size==1)|(x2.size==1)|(x1.shape==x2.shape), "I can't deal with the shapes of those arrays, Dave")


	if (x1.size == 1):
		x1 = x1 * np.ones(x2.shape);

	if (x2.size == 1):
		x2 = x2 * np.ones(x1.shape);

	opshape = x1.shape

	x1 = x1.flatten()
	x2 = x2.flatten()

	v = np.zeros(x1.size)

	#pdb.set_trace()
	signs = np.sign(x1 - x2);
	I = signs == -1;
	swap = x1[I]
	x1[I] = x2[I];
	x2[I] = swap;

	# Case 1: arguments of different signs, no problems with loss of accuracy
	I1 = (x1*x2)<0;
	# Case 2: x1 = x2
	I2 = x1 == x2;
	# Case 3: Both arguments are non-negative
	I3 = (x1 > 0) & ~I1 & ~I2;
	# Case 4: Both arguments are non-positive
	I4 = ~I1 & ~I2 & ~I3;

	#warnState = warning('query', 'MATLAB:log:logOfZero');
	#warning('off', 'MATLAB:log:logOfZero');
	#TODO in np

	v[I1] = np.log( erf(x1[I1]) - erf(x2[I1]) );
	v[I2] = -np.inf;
	v[I3] = np.log(erfcx(  x2[I3]) - erfcx(x1[I3]) * np.exp(x2[I3]**2 - x1[I3]**2)) - x2[I3]**2
	v[I4] = np.log(erfcx(  -x1[I4]) - erfcx(-x2[I4]) * np.exp(x1[I4]**2 - x2[I4]**2)) - x1[I4]**2

	#warning(warnState.state, 'MATLAB:log:logOfZero');
	#TODO

	return(v.reshape(opshape),signs.reshape(opshape))
	#wtf?if nargout < 2,
	 # v(I) = v(I) + pi*1i;
	#end

#function [h, dh_dD_i, dh_dD_j, dh_dsigma] = simComputeH(t1, t2, D_i, D_j, delta_i, delta_j, sigma)
def simComputeH(t1, t2, D_i, D_j, delta_i, delta_j, sigma, compute_gradDdecay1=False, compute_gradDdecay2=False, compute_gradsigma=False):
	"""
	(h, gradDdecay1, gradDdecay2, gradL) = simComputeH(t1, t2, D_i, D_j, delta_i, delta_j, sigma, *compute_flags)

	SIMCOMPUTEH Helper function for comptuing part of the SIM kernel.
	FORMAT
	DESC computes a portion of the SIM kernel.
	ARG t1 : first time input (number of time points x 1).
	ARG t2 : second time input (number of time points x 1).
	ARG decay1 : Decay rate for first system.
	ARG decay2 : Decay rate for second system.
	ARG l : length scale of latent process.
	RETURN h : result of this subcomponent of the kernel for the given values.
	
	DESC computes a portion of the SIM kernel and gradients with
	respect to various parameters.
	ARG t1 : first time input (number of time points x 1).
	ARG t2 : second time input (number of time points x 1).
	ARG decay1 : Decay rate for first system.
	ARG decay2 : Decay rate for second system.
	ARG l : length scale of latent process.
	RETURN h : result of this subcomponent of the kernel for the given values.
	RETURN grad_D_decay1 : gradient of H with respect to DECAY1.
	RETURN grad_D_decay2 : gradient of H with respect to DECAY1.
	RETURN grad_L : gradient of H with respect to length scale of
	latent process.
	
	COPYRIGHT : Neil D. Lawrence, 2006
	
	MODIFICATIONS : Antti Honkela, 2007, 2008
	
	MODIFICATIONS : James Hensman, 2011
	
	SEEALSO : simKernParamInit, lnDiffErfs

	KERN
	"""

	assert(((t1.shape[1]==1) and (t2.shape[1]==1)),'Input can only have one column')

	dim1 = t1.shape[0];
	dim2 = t2.shape[0];
	t1 = t1 - delta_i;
	t2 = t2 - delta_j;

	t1Mat = t1*np.ones((1,dim2))#(:, ones(1, dim2));
	t2Mat = t2*np.ones((1,dim1))#(:, ones(1, dim1))';
	t2Mat = t2Mat.T
	diffT = (t1Mat - t2Mat);

	invSigmaDiffT = 1/sigma*diffT;
	halfSigmaD_i = 0.5*sigma*D_i;
	#ind = find(t1Mat~=0); # this never gets used? -- James
	h = np.zeros(diffT.shape);

	(lnPart1, sign1) = lnDiffErfs(halfSigmaD_i + t2Mat/sigma,  halfSigmaD_i - invSigmaDiffT);
	(lnPart2, sign2) = lnDiffErfs(halfSigmaD_i, halfSigmaD_i - t1Mat/sigma);

	h = sign1 * np.exp(halfSigmaD_i*halfSigmaD_i-D_i*diffT+lnPart1 - np.log(D_i + D_j)) - sign2 * np.exp(halfSigmaD_i*halfSigmaD_i-D_i*t1Mat-D_j*t2Mat + lnPart2 - np.log(D_i + D_j));
	sigma2 = sigma*sigma;

	dh_dD_i, dh_dD_j, dh_dsigma = None, None, None

	if compute_gradDdecay1:
		dh_dD_i = (0.5*D_i*sigma2*(D_i + D_j)-1)*h + (-diffT*sign1*np.exp(halfSigmaD_i*halfSigmaD_i-D_i*diffT+lnPart1) +t1Mat*sign2*np.exp(halfSigmaD_i*halfSigmaD_i-D_i*t1Mat - D_j*t2Mat+lnPart2)) + sigma/np.sqrt(np.pi)*(-np.exp(-diffT*diffT/sigma2)+np.exp(-t2Mat*t2Mat/sigma2-D_i*t1Mat)+np.exp(-t1Mat*t1Mat/sigma2-D_j*t2Mat)-np.exp(-(D_i*t1Mat + D_j*t2Mat)))
		dh_dD_i = np.real(dh_dD_i/(D_i+D_j));
	if compute_gradDdecay2:
		dh_dD_j = t2Mat*sign2*np.exp(halfSigmaD_i*halfSigmaD_i-(D_i*t1Mat + D_j*t2Mat)+lnPart2)-h;
		dh_dD_j = np.real(dh_dD_j/(D_i + D_j));
    
	if compute_gradsigma:
		dh_dsigma = 0.5*D_i*D_i*sigma*h+ 2/(np.sqrt(np.pi)*(D_i+D_j))*((-diffT/sigma2-D_i/2)*np.exp(-diffT*diffT/sigma2)+ (-t2Mat/sigma2+D_i/2)*np.exp(-t2Mat*t2Mat/sigma2-D_i*t1Mat) - (-t1Mat/sigma2-D_i/2)*np.exp(-t1Mat*t1Mat/sigma2-D_j*t2Mat)- D_i/2*np.exp(-(D_i*t1Mat+D_j*t2Mat)));
	
	return (h, dh_dD_i, dh_dD_j, dh_dsigma)

