function g = tensorKernGradient(kern, x, varargin)

% TENSORKERNGRADIENT Gradient of TENSOR kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% tensor product
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
% SEEALSO tensorKernParamInit, kernGradient, tensorKernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

  
% Last of varargin is covGrad.
  g = zeros(1, kern.nParams);
  startVal = 1;
  endVal = 0;
  twoXin = 0;
  if length(varargin) > 1
    twoXin = 1;
    x2 = varargin{1};
    covGrad = varargin{end};
  else
    covGrad = varargin{end};
  end
  
  for i = 1:length(kern.comp)
    tempKern = tensorKernSlash(kern, i);
    if twoXin
      tempCovGrad = covGrad.*kernCompute(tempKern, x, x2);
    else
      tempCovGrad = covGrad.*kernCompute(tempKern, x);
  end
  endVal = endVal + kern.comp{i}.nParams;
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved in the kernel.
    if ~twoXin
      g(1, startVal:endVal)  = kernGradient(kern.comp{i}, ...
                                            x(:, kern.comp{i}.index), ...
                                            tempCovGrad);
    else
      g(1, startVal:endVal) = kernGradient(kern.comp{i}, ...
                                           x(:, kern.comp{i}.index), ...
                                           x2(:, kern.comp{i}.index), ...
                                           tempCovGrad);
    end
  else
    if ~twoXin
      % all the data is involved with the kernel.
      g(1, startVal:endVal)  = kernGradient(kern.comp{i}, x, ...
                                            tempCovGrad);
    else
      g(1, startVal:endVal) = kernGradient(kern.comp{i}, x, x2, ...
                                           tempCovGrad);
    end
  end
  startVal = endVal + 1;
end
g = g*kern.paramGroups;
