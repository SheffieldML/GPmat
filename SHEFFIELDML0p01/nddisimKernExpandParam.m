function kern = nddisimKernExpandParam(kern, params)

% NDDISIMKERNEXPANDPARAM Create kernel structure from NDDISIM kernel's parameters.
%
%	Description:
%
%	KERN = NDDISIMKERNEXPANDPARAM(KERN, PARAM) returns a single input
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


%fprintf(1,'Expanding DISIM parameters, received these values:\n');
%params


kern.inverseWidth = params(1);
kern.di_variance = params(2);
kern.decay = params(3);
kern.variance = params(4);
kern.delay = params(5);
