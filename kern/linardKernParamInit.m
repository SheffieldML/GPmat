function kern = linardKernParamInit(kern)

% LINARDKERNPARAMINIT LINARD kernel parameter initialisation.
% The automatic relevance determination version of the linear
% kernel (LINARD) is the simple inner product kernel with feature
% selection applied.
%
% k(x_i, x_j) = sigma2 * x_i'*A* x_j
%
% where A is a diagonal matrix of values constrained to be between
% zero and one. These parameters are stored in the field
% 'inputScales'. There is also the parameter, sigma2, which is stored
% in the field kern.variance.
%
% SEEALSO : linKernParamInit, rbfardKernParamInit
%
% FORMAT
% DESC initialises the automatic relevance determination linear
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


% This parameters is restricted positive.
kern.variance = 1;
% These parameters are restricted to lie between 0 and 1.
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = 1 + kern.inputDimension;

kern.transforms(1).index = 1;
kern.transforms(1).type = optimiDefaultConstraint('positive');
kern.transforms(2).index = [2:kern.nParams];
kern.transforms(2).type = optimiDefaultConstraint('zeroone');

kern.isStationary = false;
