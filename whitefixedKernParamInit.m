function kern = whitefixedKernParamInit(kern)

% WHITEFIXEDKERNPARAMINIT WHITEFIXED kernel parameter initialisation.
% The white noise kernel arises from assuming independent Gaussian
% noise for each point in the function. The variance of the noise is
% given by the kern.variance parameter. The fixed white noise kernel
% is a simple variant of the white kernel that doesn't allow the noise
% parameter to be optimised. It is useful when the level of noise is
% known a priori.
% 
% This kernel is not intended to be used independently, it is provided
% so that it may be combined with other kernels in a compound
% kernel.
%
% SEEALSO : cmpndKernParamInit
%
% FORMAT
% DESC initialises the fixed parameter white noise
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Nathaniel J. King, 2006

% KERN

kern.variance = exp(-2);
kern.nParams = 0;

kern.isStationary = true;
