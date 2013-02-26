function gX = noneKernDiagGradX(kern, X)

% NONEKERNDIAGGRADX Gradient of NONE kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = NONEKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the dummy kernel function kernel matrix with respect to
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
%	NONEKERNPARAMINIT, KERNDIAGGRADX, NONEKERNGRADX


%	Copyright (c) 2008 Neil D. Lawrence


gX = zeros(size(X));
