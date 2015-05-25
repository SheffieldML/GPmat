function gT = rbfinfwhiteKernDiagGradX(kern, t)

% RBFINFWHITEKERNDIAGGRADX Gradient of RBF-WHITE kernel's (with integration
% limits between minus infinity and infinity) diagonal w.r.t. t.
% FORMAT
% DESC computes the gradient of the diagonal of the RBF-WHITE kernel matrix
% with respect to the elements of the input column vector given in t.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG t : the input data in the form of a design matrix.
% RETURN gT : the gradients of the diagonal with respect to each element
% of t. The returned matrix has the same dimensions as t.
%
% SEEALSO : rbfinfwhiteKernParamInit, kernDiagGradX, rbfinfwhiteKernGradX
%
% COPYRIGHT : David Luengo, 2009

% KERN


if size(t, 2) > 1
  error('Input can only have one column');
end

gT = zeros(size(t));
