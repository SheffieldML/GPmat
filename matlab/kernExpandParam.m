function kern = kernExpandParam(kern, params)

% KERNEXPANDPARAM Expand parameters to form a kernel structure.

% IVM

kern = feval([kern.type 'KernExpandParam'], kern, params);
