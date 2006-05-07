function g = multiKernGradient(kern, x, x2, covGrad)

% MULTIKERNGRADIENT Gradient of MULTI kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% multiple output block
% kernel's parameters. As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x : the input locations for which the gradients are being
% computed. 
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations associated with the rows of the
% kernel matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO multiKernParamInit, kernGradient, multiKernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

% Collate arguments.
dim1 = size(x, 1);
arg{1} = x;
if nargin > 3
  dim2 = size(x2, 1);
  arg{2} = x2;
else
  dim2 = dim1;
  covGrad = x2;
end

g = zeros(1, size(kern.paramGroups, 1));
startVal = 1;
endVal = 0;
for i = 1:kern.numBlocks
  endVal = endVal + kern.comp{i}.nParams;

  startOne = (i-1)*dim1 + 1;
  endOne = i*dim1;
  g(1, startVal:endVal) = multiKernGradientBlock(kern, ...
                                                 arg{:}, ...
                                                 covGrad(startOne:endOne, ...
                                                    (i-1)*dim2 + 1:i*dim2), ...
                                                 i, i);
  startVal2 = 1;
  endVal2 = 0;
  for j = 1:i-1
    endVal2 = endVal2 + kern.comp{j}.nParams;
    if ~isempty(kern.block{i}.cross{j})
      startTwo = (j-1)*dim2 + 1;
      endTwo = j*dim2;
      [g1, g2] = multiKernGradientBlock(kern, ...
                                        arg{:}, ...
                                        covGrad(startOne:endOne, ...
                                                startTwo:endTwo), ...
                                        i, j);
      g(1, startVal:endVal) = g(1, startVal:endVal) + 2*g1;
      g(1, startVal2:endVal2) = g(1, startVal2:endVal2) + 2*g2;
    end
    startVal2 = endVal2 + 1;
  end
  startVal = endVal + 1;
end
g = g*kern.paramGroups;

