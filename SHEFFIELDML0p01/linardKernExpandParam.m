function kern = linardKernExpandParam(kern, params)

% LINARDKERNEXPANDPARAM Create kernel structure from LINARD kernel's parameters.
%
%	Description:
%
%	KERN = LINARDKERNEXPANDPARAM(KERN, PARAM) returns a automatic
%	relevance determination linear kernel structure filled with the
%	parameters in the given vector. This is used as a helper function to
%	enable parameters to be optimised in, for example, the NETLAB
%	optimisation functions.
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
%	LINARDKERNPARAMINIT, LINARDKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



kern.variance = params(1);
kern.inputScales = params(2:end);
