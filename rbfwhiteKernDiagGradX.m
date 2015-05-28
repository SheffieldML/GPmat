function gT = rbfwhiteKernDiagGradX(kern, t)

% RBFWHITEKERNDIAGGRADX Gradient of RBF-WHITE kernel's diagonal w.r.t. t.
% FORMAT
% DESC computes the gradient of the diagonal of the RBF-WHITE kernel matrix
% with respect to the elements of the input column vector given in t.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG t : the input data in the form of a design matrix.
% RETURN gT : the gradients of the diagonal with respect to each element
% of t. The returned matrix has the same dimensions as t.
%
% SEEALSO : rbfwhiteKernParamInit, kernDiagGradX, rbfwhiteKernGradX
%
% COPYRIGHT : David Luengo, 2009

% KERN


if size(t, 2) > 1
  error('Input can only have one column');
end

if (kern.isStationary == false)
    gT = 0.5 * kern.variance * kern.inverseWidth / pi ...
        * exp(-kern.inverseWidth * t .* t);
else
    gT = zeros(size(t));
end
