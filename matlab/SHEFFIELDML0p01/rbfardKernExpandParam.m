function kern = rbfardKernExpandParam(kern, params)

% RBFARDKERNEXPANDPARAM Create kernel structure from RBFARD kernel's parameters.
%
%	Description:
%
%	KERN = RBFARDKERNEXPANDPARAM(KERN, PARAM) returns a automatic
%	relevance determination radial basis function kernel structure
%	filled with the parameters in the given vector. This is used as a
%	helper function to enable parameters to be optimised in, for
%	example, the NETLAB optimisation functions.
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
%	RBFARDKERNPARAMINIT, RBFARDKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



kern.inverseWidth = params(1);
kern.variance = params(2);
kern.inputScales = params(3:end);
