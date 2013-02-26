function kern = dexpKernExpandParam(kern, params)

% DEXPKERNEXPANDPARAM Create a kernel structure from the double exponential
%
%	Description:
%	kernel's parameters.
%	
%
%	KERN = DEXPKERNEXPANDPARAM(KERN, PARAMS) returns a double
%	exponential kernel kernel structure filled with the parameters in
%	the given vector. This is used as a helper function to enable
%	parameters to be optimised in, for example, the NETLAB optimisation
%	functions.
%	 Returns:
%	  KERN - kernel structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  KERN - the kernel structure in which the parameters are to be
%	   placed.
%	  PARAMS - vector of parameters which are to be placed in the kernel
%	   structure.
%	
%
%	See also
%	DEXPKERNPARAMINIT, DEXPKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2009 David Luengo



kern.decay = params(1);
kern.variance = params(2);
