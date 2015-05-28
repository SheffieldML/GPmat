function kern = rbfwhiteKernExpandParam(kern, params)

% RBFWHITEKERNEXPANDPARAM Create kernel structure from RBF-WHITE kernel's
%
%	Description:
%	parameters.
%
%	KERN = RBFWHITEKERNEXPANDPARAM(KERN, PARAM) returns a RBF-WHITE
%	kernel structure filled with the parameters in the given vector.
%	This is used as a helper function to enable parameters to be
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
%	RBFWHITEKERNPARAMINIT, RBFWHITEKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2009 David Luengo



kern.inverseWidth = params(1);
kern.variance = params(2);
