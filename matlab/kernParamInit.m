function kern = kernParamInit(kern)

% KERNPARAMINIT Kernel parameter initialisation.

% KERN


kern = feval([kern.type 'KernParamInit'], kern);
