function kern = rbfardKernParamInit(kern)

% RBFARDKERNPARAMINIT RBFARD kernel parameter initialisation.
% The automatic relevance determination version of the radial basis
% function kernel (RBFARD) is a very smooth non-linear kernel and is a
% popular choice for generic use.
%
% k(x_i, x_j) = sigma2 * exp(-gamma/2 *(x_i - x_j)'*A*(x_i - x_j))
%
% The parameters are sigma2, the process variance (kern.variance), the
% diagonal matrix of input scales (kern.inputScales, constrained to be
% between zero and one) and gamma, the inverse width
% (kern.inverseWidth). The inverse width controls how wide the basis
% functions are, the larger gamma, the smaller the basis functions
% are.
%
% SEEALSO : rbfKernParamInit
%
% FORMAT
% DESC initialises the automatic relevance determination radial basis function
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


% This parameter is restricted positive.
kern.inverseWidth = 2/kern.inputDimension;
kern.variance = 1;
% These parameters are restricted to lie between 0 and 1.
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = 2 + kern.inputDimension;

kern.transforms(1).index = [1 2];
kern.transforms(1).type = optimiDefaultConstraint('positive');
kern.transforms(2).index = [3:kern.nParams];
kern.transforms(2).type = optimiDefaultConstraint('zeroone');

kern.isStationary = true;
