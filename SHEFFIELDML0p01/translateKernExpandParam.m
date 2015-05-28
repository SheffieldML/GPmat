function kern = translateKernExpandParam(kern, params)

% TRANSLATEKERNEXPANDPARAM Create kernel structure from TRANSLATE kernel's parameters.
%
%	Description:
%
%	KERN = TRANSLATEKERNEXPANDPARAM(KERN, PARAM) returns a input space
%	translation kernel structure filled with the parameters in the given
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
%	cmpndKernExpandParam, kernExpandParam
%	
%
%	See also
%	TRANSLATEKERNPARAMINIT, TRANSLATEKERNEXTRACTPARAM, 


%	Copyright (c) 2007 Neil D. Lawrence


endVal = length(params)-kern.inputDimension;
kern = cmpndKernExpandParam(kern, params(1, 1:endVal));
startVal = endVal + 1;
kern.centre = params(1, startVal:end);
