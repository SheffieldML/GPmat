function gX = cmpndKernGradX(kern, X, X2)

% CMPNDKERNGRADX Gradient of CMPND kernel with respect to a point x.
%
%	Description:
%
%	G = CMPNDKERNGRADX(KERN, X) computes the gradient of the compound
%	kernel with respect to the input positions.
%	 Returns:
%	  G - the returned gradients. The gradients are returned in a matrix
%	   which is numData x numInputs x numData. Where numData is the
%	   number of data points and numInputs is the number of input
%	   dimensions in X.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X - locations against which gradients are being computed.
%
%	G = CMPNDKERNGRADX(KERN, X1, X2) computes the gradident of the
%	compound kernel with respect to the input positions where both the
%	row positions and column positions are provided separately.
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
%	% SEEALSO CMPNDKERNPARAMINIT, KERNGRADX, CMPNDKERNDIAGGRADX


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



i = 1;
fhandle = str2func([kern.comp{i}.type 'KernGradX']);

if ~isempty(kern.comp{i}.index)
  % only part of the data is involved with the kernel.

  gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
  gX(:, kern.comp{i}.index, :) = fhandle(kern.comp{i}, ...
                                      X(:, kern.comp{i}.index), ...
                                      X2(:, kern.comp{i}.index));
else
  % all the data is involved with the kernel.
  gX = fhandle(kern.comp{i}, X, X2);
end
for i = 2:length(kern.comp)
  fhandle = str2func([kern.comp{i}.type 'KernGradX']);
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved with the kernel.
    gX(:, kern.comp{i}.index, :) = ...
        gX(:, kern.comp{i}.index, :) + ...
        fhandle(kern.comp{i}, ...
                X(:, kern.comp{i}.index), ...
                X2(:, kern.comp{i}.index));
  else
    % all the data is involved with the kernel.
    gX = gX + fhandle(kern.comp{i}, X, X2);
  end
end
