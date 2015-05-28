function kern = diagKernParamInit(kern)

% DIAGKERNPARAMINIT DIAG kernel parameter initialisation.
%
%	Description:
%	The diag covariance function takes a one dimensional input and outputs a diagonal noise that is provided by an exponentiated and scaled version of the input.
%	
%	k(x_i, x_j) = delta_ij sigma2 exp(x_i)
%	
%	The only parameter is sigma2, the process variance (kern.variance).
%	
%	
%
%	KERN = DIAGKERNPARAMINIT(KERN) initialises the diagonal noise
%	covariance function kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	WHITEKERNPARAMINIT, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2011 Neil D. Lawrence


kern.variance = exp(-2);
kern.nParams = 1;

kern.transforms.index = 1;
kern.transforms.type = optimiDefaultConstraint('positive');
kern.trans = optimiDefaultConstraint('positive');


kern.isStationary = false;
