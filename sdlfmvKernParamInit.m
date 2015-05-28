function  kern = sdlfmvKernParamInit(kern)

% SDLFMVKERNPARAMINIT SDLFMV kernel initialization
% FORMAT
% DESC
% Initializes the switching dynamical latent force model structure for
% velocities with some initial parameters. The initial parameters are 
% passed through an option in kern.
% RETURN kern : the kernel structure with the default parameters placed in.
% ARG kern : the kernel structure which requires initialisation.
%
% SEEALSO : kernCreate, kernParamInit, sdlfmKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

% Actually, this works just as a wrapper function to use the same
% parameters that the SDLFM kernel.

kern = sdlfmKernParamInit(kern);
