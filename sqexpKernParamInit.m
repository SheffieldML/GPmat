function kern = sqexpKernParamInit(kern)

% SQEXPKERNPARAMINIT SQEXP kernel parameter initialisation.
% This kernel is a 'pre-packaged' compound kernel of the form
% {'rbf', 'lin', 'bias', 'white'}. Using this kernel removes
% the overhead of mutliple calls through the 'cmpnd' kernel.
% 
% SEEALSO sqexpKernParamInit
%
% FORMAT
% DESC initialises the pre-built compound squared exponential
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004

% KERN


kern.inverseWidth = 1;
kern.rbfVariance = 1;
kern.whiteVariance = 1; 
kern.biasVariance = 1;
kern.nParams = 4;

kern.transforms(1).index = [1 2 3 4];
kern.transforms(1).type = optimiDefaultConstraint('positive');

kern.isStationary = false;
