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
%
% COPYRIGHT : Pei Gao, 2007

% MODIFICATIONS : Mauricio Alvarez, 2008

% KERN

if iscell(x)
  dim = 0;
  for i = 1:length(x)
    dim = dim + size(x{i},1);
  end
  k = zeros(dim, 1);
  startVal = 1;
  endVal = size(x{1},1);
  for i = 1:length(kern.comp)
      if ~isempty(x{i})
          k(startVal:endVal) = kernDiagCompute(kern.comp{i}, x{i});
      end
      startVal = endVal + 1;
      if i+1 <= length(kern.comp)
          endVal = endVal + size(x{i+1},1);
      end
  end
else
  k = zeros(size(x, 1)*kern.numBlocks, 1);
  endVal = size(x, 1);
  startVal = 1;
  for i = 1:length(kern.comp)
    k(startVal:endVal, 1)  = kernDiagCompute(kern.comp{i}, x);
    startVal = endVal + 1;
    endVal = endVal + size(x, 1);
  end
end

