function kern = matern32KernExpandParam(kern, params)

% MATERN32KERNEXPANDPARAM Create kernel structure from MATERN32 kernel's parameters.
%
%	Description:
%
%	KERN = MATERN32KERNEXPANDPARAM(KERN, PARAM) returns a matern kernel
%	with nu=3/2 kernel structure filled with the parameters in the given
%	vector. This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
%	 Returns:
%	  KERN - kernel structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  KERN - the kernel structure in which the parameters are to be
%	   placed.
%	  PARAM - vector of parameters which are to be placed in the kernel
%	   structure.
%	
%
%	See also
%	MATERN32KERNPARAMINIT, MATERN32KERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2006 Neil D. Lawrence


kern.lengthScale = params(1);
kern.variance = params(2);
