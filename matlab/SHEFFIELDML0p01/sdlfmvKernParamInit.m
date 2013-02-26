function  kern = sdlfmvKernParamInit(kern)

% SDLFMVKERNPARAMINIT SDLFMV kernel initialization
%
%	Description:
%
%	KERN = SDLFMVKERNPARAMINIT(KERN)  Initializes the switching
%	dynamical latent force model structure for velocities with some
%	initial parameters. The initial parameters are passed through an
%	option in kern.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	KERNCREATE, KERNPARAMINIT, SDLFMKERNPARAMINIT


%	Copyright (c) 2010 Mauricio A. Alvarez


% Actually, this works just as a wrapper function to use the same
% parameters that the SDLFM kernel.

kern = sdlfmKernParamInit(kern);
