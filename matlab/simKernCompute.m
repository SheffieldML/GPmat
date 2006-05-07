function k = simKernCompute(kern, t, t2)

% SIMKERNCOMPUTE Compute the SIM kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the single input motif
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : the input matrix associated with the rows of the kernel.
% ARG t2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the single input motif
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : simKernParamInit, kernCompute, kernCreate, simKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

if nargin < 3
  t2 = t;
end
if size(t, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

sigma = sqrt(2/kern.inverseWidth);

h = simComputeH(t, t2, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
if nargin < 3
  k = h + h';
else
  h2 = simComputeH(t2, t, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
  k = h + h2';
end
k = 0.5*k*sqrt(pi)*sigma;
k = kern.variance*k;
