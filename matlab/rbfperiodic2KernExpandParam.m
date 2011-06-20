function kern = rbfperiodic2KernExpandParam(kern, params)

% RBFPERIODIC2KERNEXPANDPARAM Create kernel structure from RBFPERIODIC2 kernel's parameters.
%
%	Description:
%
%	KERN = RBFPERIODIC2KERNEXPANDPARAM(KERN, PARAM) returns a RBF derived
%	periodic kernel structure filled with the parameters in the given
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
%	RBFPERIODIC2KERNPARAMINIT, RBFPERIODI2CKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2007 Neil D. Lawrence
%	MODIFICATIONS : Andreas C. Damianou,  Michalis K. Titsias, 2011


kern.inverseWidth = params(1);
kern.variance = params(2);
kern.factor = params(3);
kern.period = 2*pi/kern.factor;
