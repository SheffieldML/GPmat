function params = kernExtractParam(kern)

% KERNEXTRACTPARAM Extract parameters from kernel structure.

% IVM

params = feval([kern.type 'KernExtractParam'], kern);
