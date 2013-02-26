function kern = sdrbfKernExpandParam(kern, params)

% SDRBFKERNEXPANDPARAM Pass parameters from params to SDRBF kernel
%
%	Description:
%
%	KERN = SDRBFKERNEXPANDPARAM(KERN, PARAM) returns a switching
%	dynamical RBF kernel structure filled with the parameters in the
%	given vector. This is used as a helper function to enable parameters
%	to be optimised in, for example, the NETLAB optimisation functions.
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
%	SDRBFKERNPARAMINIT, SDRBFKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2010 Mauricio A. Alvarez


kern.inverseWidth = reshape(params(kern.inverseWidthIndx), kern.nlfPerInt, kern.nIntervals);
kern.switchingTimes = params(kern.switchingTimesIndx);
