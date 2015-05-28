function kern = invcmpndKernExpandParam(kern, params)

% INVCMPNDKERNEXPANDPARAM Create kernel structure from INVCMPND kernel's parameters.
%
%	Description:
%
%	KERN = INVCMPNDKERNEXPANDPARAM(KERN, PARAM) returns a inv-compound
%	kernel structure filled with the parameters in the given vector.
%	This is used as a helper function to enable parameters to be
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
%	INVCMPNDKERNPARAMINIT, INVCMPNDKERNEXTRACTPARAM, CMPNDKERNEXPANDPARAM


%	Copyright (c) 2012 Andreas C. Damianou



kern = cmpndKernExpandParam(kern, params);