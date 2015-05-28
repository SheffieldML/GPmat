function kern = matern52KernExpandParam(kern, params)

% MATERN52KERNEXPANDPARAM Create kernel structure from MATERN52 kernel's parameters.
%
%	Description:
%
%	KERN = MATERN52KERNEXPANDPARAM(KERN, PARAM) returns a matern kernel
%	with nu=5/2 kernel structure filled with the parameters in the given
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
%	MATERN52KERNPARAMINIT, MATERN52KERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2006 Neil D. Lawrence


kern.lengthScale = params(1);
kern.variance = params(2);
