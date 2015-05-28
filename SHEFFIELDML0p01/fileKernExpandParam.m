function kern = fileKernExpandParam(kern, params)

% FILEKERNEXPANDPARAM Create kernel structure from FILE kernel's parameters.
%
%	Description:
%
%	KERN = FILEKERNEXPANDPARAM(KERN, PARAM) returns a stored file kernel
%	structure filled with the parameters in the given vector. This is
%	used as a helper function to enable parameters to be optimised in,
%	for example, the NETLAB optimisation functions.
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
%	FILEKERNPARAMINIT, FILEKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2005, 2006 Neil D. Lawrence



kern.variance = params(1);
