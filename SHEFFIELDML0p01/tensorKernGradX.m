function gX = tensorKernGradX(kern, X, X2)

% TENSORKERNGRADX Gradient of TENSOR kernel with respect to a point x.
%
%	Description:
%
%	G = TENSORKERNGRADX(KERN, X) computes the gradient of the tensor
%	product kernel with respect to the input positions.
%	 Returns:
%	  G - the returned gradients. The gradients are returned in a matrix
%	   which is numData x numInputs x numData. Where numData is the
%	   number of data points and numInputs is the number of input
%	   dimensions in X.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X - locations against which gradients are being computed.
%
%	G = TENSORKERNGRADX(KERN, X1, X2) computes the gradident of the
%	tensor product kernel with respect to the input positions where both
%	the row positions and column positions are provided separately.
%	 Returns:
%	  G - the returned gradients. The gradients are returned in a matrix
%	   which is numData2 x numInputs x numData1. Where numData1 is the
%	   number of data points in X1, numData2 is the number of data points
%	   in X2 and numInputs is the number of input dimensions in X.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X1 - row locations against which gradients are being computed.
%	  X2 - column locations against which gradients are being computed.
%	
%
%	See also
%	% SEEALSO TENSORKERNPARAMINIT, KERNGRADX, TENSORKERNDIAGGRADX


%	Copyright (c) 2006 Neil D. Lawrence



if nargin > 2
  twoXin = 1;
else
  twoXin = 0;
end
i = 1;

% Create a kernel with component i missing.
tempKern = tensorKernSlash(kern, i);
% Compute kernel matrix from that kernel.
Kslash = kernCompute(tempKern, X, X2)';


if ~isempty(kern.comp{i}.index)
  % only part of the data is involved with the kernel.

  gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
  gX(:, kern.comp{i}.index, :) = kernGradX(kern.comp{i}, ...
                                      X(:, kern.comp{i}.index), ...
                                      X2(:, kern.comp{i}.index));
  gX(:, kern.comp{i}.index, :) = ...
      gX(:, kern.comp{i}.index, :) ...
      .*repmat(reshape(Kslash, [size(X2, 1) 1 size(X, 1)]), ...
               [1 length(kern.comp{i}.index) 1]);
else
  % all the data is involved with the kernel.
  gX = kernGradX(kern.comp{i}, X, X2);

  gX = gX.*repmat(reshape(Kslash, ...
                          [size(X2, 1) 1 size(X, 1)]), ...
                  [1 size(X2, 2) 1]);
  
end
for i = 2:length(kern.comp)
  % Create a kernel with component i missing.
  tempKern = tensorKernSlash(kern, i);
  % Compute kernel matrix from that kernel.
  Kslash = kernCompute(tempKern, X, X2)';
  
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved with the kernel.
    gX(:, kern.comp{i}.index, :) = ...
        gX(:, kern.comp{i}.index, :) + ...
        kernGradX(kern.comp{i}, ...
                  X(:, kern.comp{i}.index), ...
                  X2(:, kern.comp{i}.index))   ...
        .*repmat(reshape(Kslash, [size(X2, 1) 1 size(X, 1)]), ...
                 [1 length(kern.comp{i}.index) 1]);
  else
    % all the data is involved with the kernel.
    gX = gX + kernGradX(kern.comp{i}, X, X2).* ...
         repmat(reshape(Kslash, ...
                        [size(X2, 1) 1 size(X, 1)]), ...
                [1 size(X2, 2) 1]);
  end
end
