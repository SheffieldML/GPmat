function kern = linard2KernExpandParam(kern, params)

% LINARD2KERNEXPANDPARAM Create kernel structure from LINARD2 kernel's parameters.
%
%	Description:
%
%	KERN = LINARD2KERNEXPANDPARAM(KERN, PARAM) returns a automatic
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
%
%	See also
%	LINARD2KERNPARAMINIT, LINARD2KERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias


kern.inputScales = params(1:end);
