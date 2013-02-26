function kern = lfmaKernParamInit(kern)

% LFMAKERNPARAMINIT LFMA kernel parameter initialisation.
%
%	Description:
%	The kernel is designed to interoperate with the multiple output block
%	kernel so that f(t) can be inferred given several different
%	instantiations of x(t).
%	
%	The parameters (m, c, delta and k) are constrained positive.
%	
%
%	KERN = LFMAKERNPARAMINIT(KERN) initialises the latent force model
%	kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	KERNCREATE, KERNPARAMINIT, LFMKERNCOMPUTE


%	Copyright (c) 2010 Mauricio A. Alvarez


% A wrapper that calls lfmKernParamInit to initialize the parameters of the
% model.

kern = lfmKernParamInit(kern);
