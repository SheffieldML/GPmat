function gX = cmpndKernDiagGradX(kern, X)

% CMPNDKERNDIAGGRADX Gradient of CMPND kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = CMPNDKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the compound kernel matrix with respect to the elements
%	of the design matrix given in X.
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
%	CMPNDKERNPARAMINIT, KERNDIAGGRADX, CMPNDKERNGRADX


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



i = 1;
if ~isempty(kern.comp{i}.index)
  % only part of the data is involved with the kernel.
  gX = zeros(size(X));
  gX(:, kern.comp{i}.index) = kernDiagGradX(kern.comp{i}, ...
                                        X(:, kern.comp{i}.index));
else
  % all the data is involved with the kernel.
  gX = kernDiagGradX(kern.comp{i}, X);
end
for i = 2:length(kern.comp)
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved with the kernel.
    gX(:, kern.comp{i}.index) = ...
        gX(:, kern.comp{i}.index) + ...
        kernDiagGradX(kern.comp{i}, ...
                  X(:, kern.comp{i}.index));
  else
    % all the data is involved with the kernel.
    gX = gX + kernDiagGradX(kern.comp{i}, X);
  end
end
