function gX = tensorKernDiagGradX(kern, X)

% TENSORKERNDIAGGRADX Gradient of TENSOR kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = TENSORKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the tensor product kernel matrix with respect to the
%	elements of the design matrix given in X.
%	 Returns:
%	  GX - the gradients of the diagonal with respect to each element of
%	   X. The returned matrix has the same dimensions as X.
%	 Arguments:
%	  KERN - the kernel structure for which gradients are being
%	   computed.
%	  X - the input data in the form of a design matrix.
%	
%
%	See also
%	TENSORKERNPARAMINIT, KERNDIAGGRADX, TENSORKERNGRADX


%	Copyright (c) 2006 Neil D. Lawrence



i = 1;
% Create a kernel with component i missing.
tempKern = tensorKernSlash(kern, i);
% Compute kernel matrix from that kernel.
diagKslash = kernDiagCompute(tempKern, X);

if ~isempty(kern.comp{i}.index)
  
  % only part of the data is involved with the kernel.
  gX = zeros(size(X));
  gX(:, kern.comp{i}.index) = kernDiagGradX(kern.comp{i}, ...
                                        X(:, kern.comp{i}.index)) ...
      .*repmat(diagKslash, 1, length(kern.comp{i}.index));
else
  % all the data is involved with the kernel.
  gX = kernDiagGradX(kern.comp{i}, X)...
       .*repmat(diagKslash, [1 size(X, 2)]);
end
for i = 2:length(kern.comp)
  % Create a kernel with component i missing.
  tempKern = tensorKernSlash(kern, i);
  % Compute kernel matrix from that kernel.
  diagKslash = kernDiagCompute(tempKern, X);

  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved with the kernel.
    gX(:, kern.comp{i}.index) = ...
        gX(:, kern.comp{i}.index) + ...
        kernDiagGradX(kern.comp{i}, ...
                  X(:, kern.comp{i}.index)) ...
        .*repmat(diagKslash, 1, length(kern.comp{i}.index));
  else
    % all the data is involved with the kernel.
    gX = gX + kernDiagGradX(kern.comp{i}, X)...
         .*repmat(diagKslash, [1 size(X, 2)]);
  end
end
