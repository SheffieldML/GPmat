function kern = noneKernParamInit(kern)

% NONEKERNPARAMINIT NONE kernel parameter initialisation.  
% This kernel is a dummy kernel which has no parameters and returns
% zeros. It is used in multi output GP latent function setups, where for
% independent output kernels there is no associated latent function kernel.
%
% FORMAT
% DESC initialises the dummy kernel function
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2008

% KERN


kern.nParams = 0;

kern.isStationary = true;
