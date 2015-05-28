function gX = linardKernDiagGradX(kern, X)

% LINARDKERNDIAGGRADX Gradient of LINARD kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = LINARDKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the automatic relevance determination linear kernel
%	matrix with respect to the elements of the design matrix given in X.
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
%	LINARDKERNPARAMINIT, KERNDIAGGRADX, LINARDKERNGRADX


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



gX = 2*kern.variance*X.*repmat(kern.inputScales, [size(X, 1), 1]);
