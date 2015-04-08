function kern = irbfKernParamInit(kern)

% IRBFKERNPARAMINIT IRBF kernel parameter initialisation.
% The integral of the radial basis function kernel (IRBF) is a one
% dimensional kernel that expresses the integral of the RBF in one
% dimension.
%
% k(x_i, x_j) = sigma2 * exp(-gamma/2 *(x_i - x_j)'*(x_i - x_j))
%
% The parameters are sigma2, the process variance (kern.variance)
% and gamma, the inverse width (kern.inverseWidth). The inverse
% width controls how wide the basis functions are, the larger
% gamma, the smaller the basis functions are.
%
% SEEALSO : rbfKernParamInit
%
% FORMAT
% DESC initialises the integral of the RBF
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2009

% KERN

