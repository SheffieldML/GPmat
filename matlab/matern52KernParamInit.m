function kern = matern52KernParamInit(kern)

% MATERN52KERNPARAMINIT MATERN52 kernel parameter initialisation.
% The Matern class of kernels is a wide class with different
% degrees of freedom parameters. This is the specific case where nu
% = 5/2.
%
% Given 
% r = sqrt((x_i - x_j)'*(x_i - x_j))
% 
% we have
% k(x_i, x_j) = sigma2*(1+sqrt(5)*r/l+5r^2/(3l^2))*exp(-sqrt(5)*r/l)
%
% The parameters are sigma2, the process variance (kern.variance),
% and l, the length scale (kern.lengthScale).
%
% FORMAT
% DESC initialises the matern kernel with nu=5/2
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

kern.variance = 1;
kern.lengthScale = 1;
kern.nParams = 2;

kern.transforms.index = [1 2];
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = true;
