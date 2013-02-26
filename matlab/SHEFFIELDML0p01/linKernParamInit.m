function kern = linKernParamInit(kern)

% LINKERNPARAMINIT LIN kernel parameter initialisation.
%
%	Description:
%	The linear kernel (LIN) is the simple inner product
%	kernel. Sampling from this kernel produces linear functions.
%	
%	k(x_i, x_j) = sigma2 * x_i'*x_j
%	
%	There is one parameter, sigma2, which is stored in the field
%	kern.variance.
%	
%	
%
%	KERN = LINKERNPARAMINIT(KERN) initialises the linear kernel
%	structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	LINARDKERNPARAMINIT, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



kern.variance = 1;
kern.nParams = 1;

kern.transforms.index = 1;
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = false;
