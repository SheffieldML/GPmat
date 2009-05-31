function [k, sk] = simKernCompute(kern, t, t2)

% SIMKERNCOMPUTE Compute the SIM kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the single input motif
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : the input matrix associated with the rows of the kernel.
% ARG t2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
% RETURN sk : unscaled kernel matrix (i.e. only 0.5 times the sum of h's
% part).
%
% FORMAT
% DESC computes the kernel matrix for the single input motif
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
% RETURN sk : unscaled kernel matrix (i.e. only 0.5 times the sum of h's
% part).
%
% SEEALSO : simKernParamInit, kernCompute, kernCreate, simKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : David Luengo, 2009

% KERN

if nargin < 3
  t2 = t;
end
if size(t, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

sigma = sqrt(2/kern.inverseWidth);

if (kern.isStationary == false)
    h = simComputeH(t, t2, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
else
    h = simComputeHStat(t, t2, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
end
if nargin < 3
  sk = 0.5 * (h + h');
else
  if (kern.isStationary == false)
      h2 = simComputeH(t2, t, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
  else
      h2 = simComputeHStat(t2, t, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
  end
  sk = 0.5 * (h + h2');
end
if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
    k = sk;
else
    k = sqrt(pi)*sigma*sk;
end
k = kern.variance*k;
