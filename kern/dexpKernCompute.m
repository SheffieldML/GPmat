function [K, sK, n1] = dexpKernCompute(kern, x1, x2)

% DEXPKERNCOMPUTE Compute the double exponential kernel,
%
% k(x_i, x_j) = 0.5 * sigma2 * theta * exp(-theta*abs(x_i - x_j)),
%
% given the parameters (theta and sigma), and t.
%
% FORMAT
% DESC computes the kernel parameters for the double exponential kernel
% given the input matrices associated with rows and columns.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG x1 : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN K : the kernel matrix computed at the given points.
% RETURN sK: normalised kernel matrix (i.e. variance set to 1).
% RETURN n1: L1 distance between each row of x1 and x2.
%
% FORMAT
% DESC computes the kernel matrix for the double exponential kernel
% given a matrix of inputs.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG x1 : the input matrix associated with the rows and the columns.
% RETURN K : the kernel matrix computed at the given points.
% RETURN sK: normalised kernel matrix (i.e. no multiplication by the decay
% or variance of the exponential).
% RETURN n1: L1 distance between each row of x1 and x2.
%
% SEEALSO : ouKernParamInit, kernCompute, kernCreate, ouKernDiagCompute,
%
% COPYRIGHT : David Luengo, 2009
%
% COPYRIGHT : Neil D. Lawrence, 2009

% KERN


if nargin < 3
  n1 = sqrt(dist2(x1, x1));
else
  n1 = sqrt(dist2(x1, x2));
end

sK = 0.5 * exp(-kern.decay * n1);
K = sK * kern.decay * kern.variance;
