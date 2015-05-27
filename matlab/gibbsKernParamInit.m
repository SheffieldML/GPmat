function kern = gibbsKernParamInit(kern)

% GIBBSKERNPARAMINIT GIBBS kernel parameter initialisation.
% This is the non stationary kernel proposed by Mark Gibbs in his 1997
% thesis. It is similar to an RBF but has a length scale that varies
% with input location. This leads to an additional term in front of
% the kernel.
%
% Given 
% r = sqrt((x_i - x_j)'*(x_i - x_j))
% 
% we have
% k(x_i, x_j) = sigma2*Z*exp(-r^2/(l(x)*l(x) + l(x')*l(x')))
%
% where
% Z = sqrt(2*l(x)*l(x')/(l(x)*l(x) + l(x')*l(x'))
%
% The parameters are sigma2, the process variance (kern.variance),
% and the parameters of l(x) which is a function that can be specified by the user, by default an MLP is used.
%
% SEEALSO : mlpCreate, gibbsKernSetLengthScaleFunc
%
% FORMAT
% DESC initialises Mark Gibbs's non-stationary
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

kern.variance = 1;
options = mlpOptions(5);
kern.lengthScaleFunc = modelCreate('mlp', kern.inputDimension, 1, options);
kern.lengthScaleTransform = optimiDefaultConstraint('positive');

kern.nParams = 1+kern.lengthScaleFunc.numParams;

kern.transforms.index = kern.nParams;
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = false;
