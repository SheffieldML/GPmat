function gX = linKernDiagGradX(kern, X)

% LINKERNDIAGGRADX Gradient of LIN kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = LINKERNDIAGGRADX(KERN, X) computes the gradient of the diagonal
%	of the linear kernel matrix with respect to the elements of the
%	design matrix given in X.
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
%	LINKERNPARAMINIT, KERNDIAGGRADX, LINKERNGRADX


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



gX = 2*X*kern.variance;
