function [K, sK] = indexardKernCompute(kern, x, x2)

% INDEXARDKERNCOMPUTE Compute the INDEXARD kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the index ard based covariance function
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the index ard based covariance function
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : indexardKernParamInit, kernCompute, kernCreate, indexardKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN
  if size(x, 2)>1
    error('Index kernel requires 1-dimensional input.')
  end

  if nargin<3
    x2 = x;
  end
  K = zeros(size(x, 1), size(x2, 1));
  for i = 1:size(x, 1)
    for j = 1:size(x2, 1)
      if round(x(i)) == round(x2(j))
        ind = find(round(x(i))==kern.indices);
        if isempty(ind)
          error('Unknown index in input');
        end
        K(i, j) = kern.indexScales(ind);
      end
    end
  end
end
