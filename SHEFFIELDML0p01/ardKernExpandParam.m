function kern = ardKernExpandParam(kern, params)

% ARDKERNEXPANDPARAM Create kernel structure from ARD kernel's parameters.
%
%	Description:
%
%	KERN = ARDKERNEXPANDPARAM(KERN, PARAM) returns a pre-built RBF and
%	linear ARD kernel structure filled with the parameters in the given
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
%	ARDKERNPARAMINIT, ARDKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2004 Neil D. Lawrence



kern.inverseWidth = params(1);
kern.rbfVariance = params(2);
kern.biasVariance = params(3);
kern.whiteVariance = params(4);
kern.linearVariance = params(5);

kern.inputScales = params(6:end);
