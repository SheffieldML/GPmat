function gX = matern52KernDiagGradX(kern, X)

% MATERN52KERNDIAGGRADX Gradient of MATERN52 kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = MATERN52KERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the matern kernel with nu=5/2 kernel matrix with respect
%	to the elements of the design matrix given in X.
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
%	MATERN52KERNPARAMINIT, KERNDIAGGRADX, MATERN52KERNGRADX


%	Copyright (c) 2006 Neil D. Lawrence


gX = zeros(size(X));
