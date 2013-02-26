function kern = noneKernParamInit(kern)

% NONEKERNPARAMINIT NONE kernel parameter initialisation.
%
%	Description:
%	This kernel is a dummy kernel which has no parameters and returns
%	zeros. It is used in multi output GP latent function setups, where for
%	independent output kernels there is no associated latent function kernel.
%	
%
%	KERN = NONEKERNPARAMINIT(KERN) initialises the dummy kernel function
%	kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2008 Neil D. Lawrence



kern.nParams = 0;

kern.isStationary = true;