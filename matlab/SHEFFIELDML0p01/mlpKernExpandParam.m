function kern = mlpKernExpandParam(kern, params)

% MLPKERNEXPANDPARAM Create kernel structure from MLP kernel's parameters.
%
%	Description:
%
%	KERN = MLPKERNEXPANDPARAM(KERN, PARAM) returns a multi-layer
%	perceptron kernel structure filled with the parameters in the given
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
%	MLPKERNPARAMINIT, MLPKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence




kern.weightVariance = params(1);
kern.biasVariance = params(2);
kern.variance = params(3);
