function kern = ratquadKernExpandParam(kern, params)

% RATQUADKERNEXPANDPARAM Create kernel structure from RATQUAD kernel's parameters.
%
%	Description:
%
%	KERN = RATQUADKERNEXPANDPARAM(KERN, PARAM) returns a rational
%	quadratic kernel structure filled with the parameters in the given
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
%	RATQUADKERNPARAMINIT, RATQUADKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2006 Neil D. Lawrence


kern.alpha = params(1);
kern.lengthScale = params(2);
kern.variance = params(3);
