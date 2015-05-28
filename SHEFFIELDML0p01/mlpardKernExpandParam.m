function kern = mlpardKernExpandParam(kern, params)

% MLPARDKERNEXPANDPARAM Create kernel structure from MLPARD kernel's parameters.
%
%	Description:
%
%	KERN = MLPARDKERNEXPANDPARAM(KERN, PARAM) returns a automatic
%	relevance determination multi-layer perceptron kernel structure
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
%	MLPARDKERNPARAMINIT, MLPARDKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



kern.weightVariance = params(1);
kern.biasVariance = params(2);
kern.variance = params(3);
kern.inputScales = params(4:end);
