function kern = gibbsKernExpandParam(kern, params)

% GIBBSKERNEXPANDPARAM Create kernel structure from GIBBS kernel's parameters.
%
%	Description:
%
%	KERN = GIBBSKERNEXPANDPARAM(KERN, PARAM) returns a Mark Gibbs's
%	non-stationary kernel structure filled with the parameters in the
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
%	GIBBSKERNPARAMINIT, GIBBSKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2006 Neil D. Lawrence


kern.lengthScaleFunc = modelExpandParam(kern.lengthScaleFunc, params(1:end-1));
kern.variance = params(end);
