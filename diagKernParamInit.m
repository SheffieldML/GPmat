function kern = diagKernParamInit(kern)

% DIAGKERNPARAMINIT DIAG kernel parameter initialisation.
% The diag covariance function takes a one dimensional input and outputs a diagonal noise that is provided by an exponentiated and scaled version of the input.
%
% k(x_i, x_j) = delta_ij sigma2 exp(x_i)
%
% The only parameter is sigma2, the process variance (kern.variance).
%
% SEEALSO : whiteKernParamInit
%
% FORMAT
% DESC initialises the diagonal noise covariance function
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN

kern.variance = exp(-2);
kern.nParams = 1;

kern.transforms.index = 1;
kern.transforms.type = optimiDefaultConstraint('positive');
kern.trans = optimiDefaultConstraint('positive');


kern.isStationary = false;
