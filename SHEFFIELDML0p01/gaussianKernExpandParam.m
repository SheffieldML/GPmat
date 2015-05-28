function kern = gaussianKernExpandParam(kern, params)

% GAUSSIANKERNEXPANDPARAM Create kernel structure from gaussian kernel's parameters.
%
%	Description:
%
%	KERN = GAUSSIANKERNEXPANDPARAM(KERN, PARAM) returns a gaussian
%	kernel structure filled with the parameters in the given vector.
%	This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
%	 Returns:
%	  KERN - kernel structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  KERN - the kernel structure in which the parameters are to be
%	  PARAM - vector of parameters which are to be placed in the kernel
%	   structure.
%	
%	
%
%	See also
%	GAUSSIANKERNPARAMINIT, GAUSSIANKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2008 Mauricio Alvarez and Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009

  
kern.sigma2Latent = params(end);
kern.precisionU =  params(1:end-1)';