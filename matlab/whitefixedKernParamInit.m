function kern = whitefixedKernParamInit(kern)

% WHITEFIXEDKERNPARAMINIT white noise kernel parameter initialisation.
% FORMAT
% DESC initialises the parameters of a fixed white kernel.
% ARG kern : the kernel structure for which the parameters will be
% initialised.
% RETURN kern : the kernel structure with the parameters
% initialised.
%
% SEEALSO : kernParamInit, whitefixedKernCompute
%
% COPYRIGHT : Nathaniel J. King, 2006

% KERN

kern.variance = exp(-2);
kern.nParams = 0;

