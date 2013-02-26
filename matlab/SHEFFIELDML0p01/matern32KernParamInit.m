function kern = matern32KernParamInit(kern)

% MATERN32KERNPARAMINIT MATERN32 kernel parameter initialisation.
%
%	Description:
%	The Matern class of kernels is a wide class with different
%	degrees of freedom parameters. This is the specific case where nu
%	= 3/2.
%	
%	Given
%	r = sqrt((x_i - x_j)'*(x_i - x_j))
%	
%	We have
%	k(x_i, x_j) = sigma2*(1+sqrt(3)*r/l)*exp(-sqrt(3)*r/l)
%	
%	The parameters are sigma2, the process variance (kern.variance),
%	and l, the length scale (kern.lengthScale).
%
%	KERN = MATERN32KERNPARAMINIT(KERN) initialises the matern kernel
%	with nu=3/2 kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2006 Neil D. Lawrence


kern.variance = 1;
kern.lengthScale = 1;
kern.nParams = 2;

kern.transforms.index = [1 2];
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = true;