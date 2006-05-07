function k = multiKernDiagCompute(kern, x)

% MULTIKERNDIAGCOMPUTE Compute diagonal of MULTI kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the multiple output block kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : multiKernParamInit, kernDiagCompute, kernCreate, multiKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

k = zeros(size(x, 1)*kern.numBlocks, 1);
endVal = size(x, 1);
startVal = 1;
for i = 1:length(kern.comp)
  k(startVal:endVal, 1)  = kernDiagCompute(kern.comp{i}, x);
  startVal = endVal + 1;
  endVal = endVal + size(x, 1);
end

