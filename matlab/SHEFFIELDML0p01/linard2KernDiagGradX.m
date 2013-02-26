function gX = linard2KernDiagGradX(kern, X)

% LINARD2KERNDIAGGRADX Gradient of LINARD2 kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = LINARD2KERNDIAGGRADX(KERN, X) computes the gradient of the
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
%
%	See also
%	LINARD2KERNPARAMINIT, KERNDIAGGRADX, LINARD2KERNGRADX


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias



gX = 2*X.*repmat(kern.inputScales, [size(X, 1), 1]);
