function kern = polyardKernParamInit(kern)

% POLYARDKERNPARAMINIT POLYARD kernel parameter initialisation.
% The automatic relevance determination version of the polynomial
% kernel is included for completeness, but its use is generally not
% recommended.
%
%  k(x_i, x_j) = sigma2*(w*x_i'*A*x_j+b)^d
%
% The kernel parameters are sigma2 (kern.variance), w
% (kern.weightVariance), b (kern.biasVariance), A
% (kern.inputScales) a diagonal matrix of input scales and d
% (kern.degree). Only gradients of the first four are provided for
% kernel optimisation, it is assumed that polynomial degree would
% be set by hand.
%
% The kernel is not recommended as it is badly behaved when the
% w*x_i'*A*x_j + b has a magnitude greater than one. 
%
% SEEALSO : polyKernParamInit
%
% FORMAT
% DESC initialises the automatic relevance determination polynomial
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


kern.weightVariance = 1;
kern.biasVariance = 1;
kern.variance = 1;
% These parameters are restricted to lie between 0 and 1.
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = 3 + kern.inputDimension;

kern.degree = 2;

kern.transforms(1).index = [1 2 3];
kern.transforms(1).type = optimiDefaultConstraint('positive');
kern.transforms(2).index = [4:kern.nParams];
kern.transforms(2).type = optimiDefaultConstraint('zeroone');

kern.isStationary = false;
