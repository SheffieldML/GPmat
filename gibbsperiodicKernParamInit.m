function kern = gibbsperiodicKernParamInit(kern)

% GIBBSPERIODICKERNPARAMINIT GIBBSPERIODIC kernel parameter initialisation.
% This kernel is a periodic kernel constructed by mapping a one
% dimensional input into a two dimensional space,
%
% u_1(x) = cos(x), u_2(x) = sin(x)
%
% A Gibbs kernel is then applied in the resulting
% two dimensional space. The resulting form of the covariance is
% then
% Given 
% r = 2*sin((x_i-x_j)/2)
% 
% we have
% k(x_i, x_j) = sigma2*Z*exp(-r^2/(l(x)*l(x) + l(x')*l(x')))
%
% where
% Z = sqrt(2*l(x)*l(x')/(l(x)*l(x) + l(x')*l(x'))
%
% The parameters are sigma2, the process variance (kern.variance),
% and the parameters of l(x) which is a periodic function that can be specified by the user.
%
% SEEALSO : mlpCreate, gibbsperiodicKernSetLengthScaleFunc%
%
% SEEALSO : gibssKernParamInit
%
% FORMAT
% DESC initialises the Gibbs-kernel derived periodic kernel structure
% with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN


kern.variance = 1;
options = rbfperiodicOptions(5);
kern.lengthScaleFunc = modelCreate('rbfperiodic', kern.inputDimension, 1, options);
kern.lengthScaleTransform = optimiDefaultConstraint('positive');

kern.nParams = 1+kern.lengthScaleFunc.numParams;

kern.transforms.index = kern.nParams;
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = false;
