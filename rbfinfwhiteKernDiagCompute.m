function k = rbfinfwhiteKernDiagCompute(kern, t)

% RBFINFWHITEKERNDIAGCOMPUTE Compute diagonal of RBF-WHITE kernel (with
% integration limits between minus infinity and infinity).
% FORMAT
% DESC computes the diagonal of the kernel matrix for the RBF-WHITE kernel
% given a column vector of inputs.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG t : input data in the form of a column vector.
% RETURN k : a vector of the same size as t containing the diagonal of the
% kernel matrix computed at the given points.
%
% SEEALSO : rbfinfwhiteKernParamInit, kernDiagCompute, kernCreate,
% rbfinfwhiteKernCompute
%
% COPYRIGHT : David Luengo, 2009

% KERN


if size(t, 2) > 1
  error('Input can only have one column');
end

k = (0.5 * kern.variance * sqrt(kern.inverseWidth/pi)) * ones(size(t));
