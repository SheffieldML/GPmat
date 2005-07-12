function kern = kernParamInit(kern)

% KERNPARAMINIT Kernel parameter initialisation.

% KERN

fhandle = str2func([kern.type 'KernParamInit']);
% By default don't transform kernel parameters.
kern.transforms = [];
kern = fhandle(kern);
