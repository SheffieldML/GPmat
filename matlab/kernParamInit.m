function kern = kernParamInit(kern)

% KERNPARAMINIT Kernel parameter initialisation.

% IVM

kern = feval([kern.type 'KernParamInit'], kern);