function kern = rbfperiodicKernExpandParam(kern, params)

% RBFPERIODICKERNEXPANDPARAM Create kernel structure from RBFPERIODIC kernel's parameters.
%
%	Description:
%
%	KERN = RBFPERIODICKERNEXPANDPARAM(KERN, PARAM) returns a RBF derived
%	periodic kernel structure filled with the parameters in the given
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
%	RBFPERIODICKERNPARAMINIT, RBFPERIODICKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2007 Neil D. Lawrence



kern.inverseWidth = params(1);
kern.variance = params(2);
