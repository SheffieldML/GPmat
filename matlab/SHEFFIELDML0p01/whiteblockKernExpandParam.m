function kern = whiteblockKernExpandParam(kern, params)

% WHITEBLOCKKERNEXPANDPARAM Fill WHITEBLOCK kernel structure with params.
%
%	Description:
%
%	KERN = WHITEBLOCKKERNEXPANDPARAM(KERN, PARAM) returns a white noise
%	block kernel structure filled with the parameters in the given
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
%	WHITEBLOCKKERNPARAMINIT


%	Copyright (c) 2010 Mauricio A. Alvarez


kern.variance = params;
