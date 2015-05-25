function g = translateKernGradient(kern, varargin)

% TRANSLATEKERNGRADIENT Gradient of TRANSLATE kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% input space translation
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
% SEEALSO translateKernParamInit, kernGradient,
% translateKernDiagGradient, cmpndKernGradient, cmpndKernGradX, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

for i = 1:length(varargin)-1
  varargin{i} = varargin{i} - repmat(kern.centre, size(varargin{i}, 1), 1);
end
g = cmpndKernGradient(kern, varargin{:});


% Gradients with respect to X.
if length(varargin) == 2
  gKX = cmpndKernGradX(kern, varargin{1}, varargin{1});
  gKX = gKX*2;
  dgKX = cmpndKernDiagGradX(kern, varargin{1});
  for i = 1:size(gKX, 1)
    gKX(i, :, i) = dgKX(i, :);
  end
else
  gKX_12 = cmpndKernGradX(kern, varargin{1}, varargin{2});
  gKX_21 = cmpndKernGradX(kern, varargin{2}, varargin{1});
end
  
gcentre = zeros(1, kern.inputDimension);
if length(varargin) == 2
  for i = 1:size(varargin{1}, 1);
    for j = 1:kern.inputDimension
      gcentre(1, j) = gcentre(1, j) - gKX(:, j, i)'*varargin{end}(:, i);
    end
  end
else
  for i = 1:size(varargin{1}, 1)
    for j = 1:kern.inputDimension
      gcentre(1, j) = gcentre(1, j) - sum(gKX_12(:, j, i).*varargin{end}(i, :)');
    end
  end
  for i = 1:size(varargin{2}, 1)
    for j = 1:kern.inputDimension
      gcentre(1, j) = gcentre(1, j) - sum(gKX_21(:, j, i).*varargin{end}(:, i));
    end
  end
end
g = [g gcentre];
