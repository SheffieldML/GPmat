function kern = lfmaKernParamInit(kern)

% LFMAKERNPARAMINIT LFMA kernel parameter initialisation. 
% The kernel is designed to interoperate with the multiple output block
% kernel so that f(t) can be inferred given several different
% instantiations of x(t).
%
% The parameters (m, c, delta and k) are constrained positive.
%
% FORMAT
% DESC initialises the latent force model kernel structure with some
% default parameters.
% RETURN kern : the kernel structure with the default parameters placed in.
% ARG kern : the kernel structure which requires initialisation.
%
% SEEALSO : kernCreate, kernParamInit, lfmKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

% A wrapper that calls lfmKernParamInit to initialize the parameters of the
% model.

kern = lfmKernParamInit(kern);
