function kern = velotransKernExpandParam(kern, params)

% VELOTRANSKERNEXPANDPARAM Create kernel structure from VELOTRANS kernel's parameters.
%
%	Description:
%
%	KERN = VELOTRANSKERNEXPANDPARAM(KERN, PARAM) returns a velocity
%	translate kernel structure filled with the parameters in the given
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
%	VELOTRANSKERNPARAMINIT, VELOTRANSKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2011 Neil D. Lawrence


endVal = length(params)-kern.inputDimension+1;
kern = cmpndKernExpandParam(kern, params(1, 1:endVal));
startVal = endVal + 1;
kern.velocity = params(1, startVal:end);
