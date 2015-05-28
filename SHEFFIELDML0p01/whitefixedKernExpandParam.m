function kern = whitefixedKernExpandParam(kern, params)

% WHITEFIXEDKERNEXPANDPARAM Create kernel structure from WHITEFIXED kernel's parameters.
%
%	Description:
%
%	KERN = WHITEFIXEDKERNEXPANDPARAM(KERN, PARAM) returns a fixed
%	parameter white noise kernel structure filled with the parameters in
%	the given vector. This is used as a helper function to enable
%	parameters to be optimised in, for example, the NETLAB optimisation
%	functions.
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
%	WHITEFIXEDKERNPARAMINIT, WHITEFIXEDKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2006 Nathaniel J. King



kern = kern;