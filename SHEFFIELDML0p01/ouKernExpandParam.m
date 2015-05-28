function kern = ouKernExpandParam(kern, params)

% OUKERNEXPANDPARAM Create kernel structure from OU kernel's parameters
%
%	Description:
%	(see ouKernCompute or ouKernParamInit for a more detailed description of
%	the OU kernel).
%
%	KERN = OUKERNEXPANDPARAM(KERN, PARAM) returns a Ornstein-Uhlenbeck
%	kernel kernel structure filled with the parameters in the given
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
%	OUKERNPARAMINIT, OUKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2009 David Luengo



kern.decay = params(1);
kern.variance = params(2);
