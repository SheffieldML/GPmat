function kern = kernExpandParam(params, kern)

% KERNEXPANDPARAM Expand parameters to form a kernel structure.

% IVM

kern = feval([kern.type 'KernExpandParam'], params,  kern);
