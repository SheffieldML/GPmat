function gX = rbfardjitKernDiagGradX(kern, X)

% RBFARDJITKERNDIAGGRADX Gradient of RBFARDJIT kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = RBFARDJITKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the automatic relevance determination radial basis
%	function kernel matrix with respect to the elements of the design
%	matrix given in X.
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
%	RBFARD2KERNPARAMINIT, KERNDIAGGRADX, RBFARD2KERNGRADX


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias



gX = zeros(size(X));