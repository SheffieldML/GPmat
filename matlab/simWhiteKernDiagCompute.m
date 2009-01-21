function k = simWhiteKernDiagCompute(kern, t)

% SIMWHITEKERNDIAGCOMPUTE Compute diagonal of SIM-WHITE kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the SIM-White (Single
% Input Motif - White) kernel given a column vector of inputs.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG t : input data in the form of a column vector.
% RETURN k : a vector of the same size as t containing the diagonal of the
% kernel matrix computed at the given points.
%
% SEEALSO : simWhiteKernParamInit, kernDiagCompute, kernCreate, simWhiteKernCompute
%
% COPYRIGHT : David Luengo, 2009

% KERN


if size(t, 2) > 1
  error('Input can only have one column');
end

c = kern.variance * (kern.sensitivity^2) / (2*kern.decay);
k = ones(size(t));
if (kern.isStationary == false)
    k = k - exp(-2*kern.decay*t);
end
k = c*k;
