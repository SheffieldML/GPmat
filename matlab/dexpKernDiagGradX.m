function gX = dexpKernDiagGradX(kern, x)

% DEXPKERNDIAGGRADX Gradient of the double exponential kernel's diagonal
% with respect to x.
%
% FORMAT
% DESC computes the gradient of the diagonal of the double exponential
% kernel matrix with respect to the elements of the inputs in x.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG x : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of x. The returned vector has the same dimensions as x.
%
% SEEALSO : dexpKernParamInit, kernDiagGradX, dexpkernGradX
%
% COPYRIGHT : David Luengo, 2009

% KERN


gX = zeros(size(x));
