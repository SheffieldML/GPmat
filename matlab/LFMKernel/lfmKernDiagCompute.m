function k = lfmKernDiagCompute(kern, t)

% LFMKERNDIAGCOMPUTE Compute diagonal of LFM kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the single input
% motif kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : lfmKernParamInit, kernDiagCompute, kernCreate, lfmKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2007

% LFM

if size(t, 2) > 1 
  error('Input can only have one column');
end

k = zeros(size(t, 1), 1);
for i = 1:size(t, 1)
  k(i) = lfmXlfmKernCompute(kern, kern, t(i), t(i));
end