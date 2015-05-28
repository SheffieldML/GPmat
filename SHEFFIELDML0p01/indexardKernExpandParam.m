function kern = indexardKernExpandParam(kern, params)

% INDEXARDKERNEXPANDPARAM Create kernel structure from INDEXARD kernel's parameters.
%
%	Description:
%
%	KERN = INDEXARDKERNEXPANDPARAM(KERN, PARAM) returns a index based
%	covariance function kernel structure filled with the parameters in
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
%	INDEXARDKERNPARAMINIT, INDEXARDKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2011 Neil D. Lawrence


  kern.indexScales = params(1:end);
end