function gX = wienerKernDiagGradX(kern, X, X2)

% WIENERKERNDIAGGRADX Gradient of WIENER kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = WIENERKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the wiener kernel matrix with respect to the elements of
%	the design matrix given in X.
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
%	WIENERKERNPARAMINIT, KERNDIAGGRADX, WIENERKERNGRADX


%	Copyright (c) 2009 Neil D. Lawrence


gX = repmat(kern.variance, [1 1 size(X, 1)]);
