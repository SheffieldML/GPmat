function kern = disimKernExpandParam(kern, params)

% DISIMKERNEXPANDPARAM Create kernel structure from DISIM kernel's parameters.
%
%	Description:
%
%	KERN = DISIMKERNEXPANDPARAM(KERN, PARAM) returns a single input
%	motif kernel structure filled with the parameters in the given
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
%
%	See also
%	DISIMKERNPARAMINIT, DISIMKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2006 Neil D. Lawrence
%	Copyright (c) 2007-2009 Antti Honkela


kern.di_decay = params(1);
kern.inverseWidth = params(2);
kern.di_variance = params(3);
kern.decay = params(4);
kern.variance = params(5);
kern.rbf_variance = params(6);
if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  kern.initialVariance = params(7);
end
