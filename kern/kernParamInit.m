function kern = kernParamInit(kern)

% KERNPARAMINIT Kernel parameter initialisation.
% FORMAT
% DESC initialises the parameters of a kernel.
% ARG kern : the kernel structure for which the parameters will be
% initialised.
% RETURN kern : the kernel structure with the parameters
% initialised.
%
% SEEALSO : kernCreate

% KERN

fhandle = str2func([kern.type 'KernParamInit']);
% By default don't transform kernel parameters.
kern.transforms = [];
kern = fhandle(kern);
