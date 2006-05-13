function K = multiKernCompute(kern, varargin)

% MULTIKERNCOMPUTE Compute the MULTI kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the multiple output block
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the multiple output block
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : multiKernParamInit, kernCompute, kernCreate, multiKernDiagCompute

% KERN

dim1 = size(varargin{1}, 1);
if length(varargin)>1
  dim2 = size(varargin{2}, 1);
else
  dim2 = dim1;
end
K = zeros(kern.numBlocks*dim1, kern.numBlocks*dim2);
for i = 1:kern.numBlocks
  startOne = (i-1)*dim1 + 1;
  endOne = i*dim1;
  startThree = (i-1)*dim2 + 1;
  endThree = i*dim2;
  K(startOne:endOne, startThree:endThree) = multiKernComputeBlock(kern, ...
                                                    varargin{:}, i, ...
                                                    i);
  for j = 1:i-1
    if ~isempty(kern.block{i}.cross{j})
      startTwo = (j-1)*dim2 + 1;
      endTwo = j*dim2;
      K(startOne:endOne, startTwo:endTwo) = multiKernComputeBlock(kern, ...
                                                        varargin{:}, ...
                                                        i, j);
      if length(varargin)<2
        K(startTwo:endTwo, startOne:endOne) = K(startOne:endOne, ...
                                                startTwo:endTwo)';
      else
        startFour = (j-1)*dim1 + 1;
        endFour = j*dim1;
        K(startFour:endFour, startThree:endThree) = ...
            multiKernComputeBlock(kern, varargin{end:-1:1}, j, i)';
      end
    end
  end
end

  