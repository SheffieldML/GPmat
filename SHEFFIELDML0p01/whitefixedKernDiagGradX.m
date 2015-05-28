function gX = whitefixedKernDiagGradX(kern, X)

% WHITEFIXEDKERNDIAGGRADX Gradient of WHITEFIXED kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = WHITEFIXEDKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the fixed parameter white noise kernel matrix with
%	respect to the elements of the design matrix given in X.
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
%	WHITEFIXEDKERNPARAMINIT, KERNDIAGGRADX, WHITEFIXEDKERNGRADX


%	Copyright (c) 2006 Nathaniel J. King



gX = whiteKernDiagGradX(kern, X);