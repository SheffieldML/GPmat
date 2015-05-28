function kern = whitehKernParamInit(kern)

% WHITEHKERNPARAMINIT WHITEH kernel parameter initialisation.
% The whiteh noise kernel arises from assuming independent Gaussian
% noise for each point in the function. The variance of the noise is
% given by the kern.variance parameter.
% 
% This kernel is not intended to be used independently, it is provided
% so that it may be combined with other kernels in a compound kernel.
%
% SEEALSO : cmpndKernParamInit
%
% FORMAT
% DESC initialises the whiteh noise
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN

kern.variance = exp(-2);
kern.nParams = 1;

kern.transforms.index = 1;
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = true;
