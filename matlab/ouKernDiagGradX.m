function gT = ouKernDiagGradX(kern, t)

% OUKERNDIAGGRADX Gradient of OU kernel's diagonal with respect to t (see
% ouKernCompute or ouKernParamInit for a more detailed description of the
% OU kernel).
% FORMAT
% DESC computes the gradient of the diagonal of the Ornstein-Uhlenbeck
% kernel matrix with respect to the elements of the column vector in t.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG t : the input data in the form of a design matrix.
% RETURN gT : the gradients of the diagonal with respect to each element
% of t. The returned vector has the same dimensions as t.
%
% SEEALSO : ouKernParamInit, kernDiagGradX, oukernGradX
%
% COPYRIGHT : David Luengo, 2009

% KERN


if (kern.isStationary == true)
    gT = zeros(size(t));
else
    gT = kern.variance*exp(-2*kern.decay*t);
end
