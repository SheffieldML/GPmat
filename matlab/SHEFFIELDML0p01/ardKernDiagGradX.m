function gX = ardKernDiagGradX(kern, X)

% ARDKERNDIAGGRADX Gradient of ARD kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = ARDKERNDIAGGRADX(KERN, X) computes the gradient of the diagonal
%	of the pre-built RBF and linear ARD kernel matrix with respect to
%	the elements of the design matrix given in X.
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
%	ARDKERNPARAMINIT, KERNDIAGGRADX, ARDKERNGRADX


%	Copyright (c) 2004 Neil D. Lawrence



gX = 2*kern.linearVariance*X.*repmat(kern.inputScales, [size(X, 1) 1]);

