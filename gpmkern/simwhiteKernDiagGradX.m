function gT = simwhiteKernDiagGradX(kern, t)

% SIMWHITEKERNDIAGGRADX Gradient of SIM-WHITE kernel's diagonal w.r.t. t.
% FORMAT
% DESC computes the gradient of the diagonal of the SIM-White (Single Input
% Motif - White) kernel matrix with respect to the elements of the input
% column vector given in t.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG t : the input data in the form of a design matrix.
% RETURN gT : the gradients of the diagonal with respect to each element
% of t. The returned matrix has the same dimensions as t.
%
% SEEALSO : simwhiteKernParamInit, kernDiagGradX, simwhitekernGradX
%
% COPYRIGHT : David Luengo, 2009

% KERN


if (kern.isStationary == true)
    gT = zeros(size(t));
else
    gT = kern.variance * (kern.sensitivity^2) * exp(-2*kern.decay*t);
end
