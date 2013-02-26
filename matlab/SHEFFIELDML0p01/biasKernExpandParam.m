function kern = biasKernExpandParam(kern, params)

% BIASKERNEXPANDPARAM Create kernel structure from BIAS kernel's parameters.
%
%	Description:
%
%	KERN = BIASKERNEXPANDPARAM(KERN, PARAM) returns a bias kernel
%	structure filled with the parameters in the given vector. This is
%	used as a helper function to enable parameters to be optimised in,
%	for example, the NETLAB optimisation functions.
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
%	BIASKERNPARAMINIT, BIASKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



kern.variance = params(1);