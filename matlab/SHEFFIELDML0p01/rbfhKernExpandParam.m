function kern = rbfhKernExpandParam(kern, params)

% RBFHKERNEXPANDPARAM Create kernel structure from RBFH kernel's parameters.
%
%	Description:
%
%	KERN = RBFHKERNEXPANDPARAM(KERN, PARAM) returns a radial basis
%	function HEAT kernel structure filled with the parameters in the
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
%	RBFHKERNPARAMINIT, RBFHKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2010 Mauricio A. Alvarez


kern.inverseWidthTime = params(1);
kern.inverseWidthSpace = params(2);
