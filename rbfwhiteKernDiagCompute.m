function k = rbfwhiteKernDiagCompute(kern, t)

% RBFWHITEKERNDIAGCOMPUTE Compute diagonal of RBF-WHITE kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the RBF-WHITE kernel
% given a column vector of inputs.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG t : input data in the form of a column vector.
% RETURN k : a vector of the same size as t containing the diagonal of the
% kernel matrix computed at the given points.
%
% SEEALSO : rbfwhiteKernParamInit, kernDiagCompute, kernCreate,
% rbfwhiteKernCompute
%
% COPYRIGHT : David Luengo, 2009

% KERN


if size(t, 2) > 1
  error('Input can only have one column');
end

if (kern.isStationary == false)
    k = kern.variance/sqrt(8*pi) * erf(kern.inverseWidth*t/sqrt(2));
else
    k = zeros(size(t));
end
