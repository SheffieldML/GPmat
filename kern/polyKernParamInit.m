function kern = polyKernParamInit(kern)

% POLYKERNPARAMINIT POLY kernel parameter initialisation.
% Included for completeness, but generally not recommended, is the
% polynomial kernel,
%
%  k(x_i, x_j) = sigma2*(w*x_i'*x_j+b)^d
%
% The kernel parameters are sigma2 (kern.variance), w
% (kern.weightVariance), b (kern.biasVariance) and d
% (kern.degree). Only gradients of the first three are provided for
% kernel optimisation, it is assumed that polynomial degree would
% be set by hand.
%
% The kernel is not recommended as it is badly behaved when the
% w*x_i'*x_j + b has a magnitude greater than one. For completeness
% there is an automatic relevance determination version of this
% kernel provided.
%
% SEEALSO : polyardKernParamInit
%
% FORMAT
% DESC initialises the polynomial
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


kern.variance = 1;
kern.weightVariance = 1;
kern.biasVariance = 1;
kern.degree = 2;
kern.nParams = 3;

kern.transforms.index = [1 2 3];
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = false;
