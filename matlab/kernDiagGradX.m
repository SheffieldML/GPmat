function k = kernDiagGradX(kern, x)

% KERNDIAGGRADX Compute the gradient of the  kernel wrt X.
% FORMAT
% DESC computes the gradient of the diagonal of the kernel matrix
% with respect to the elements of the design matrix given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each
% element of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : kernDiagGradX, kernGradX

% KERN

fhandle = str2func([kern.type 'KernDiagGradX']);
k = fhandle(kern, x);
