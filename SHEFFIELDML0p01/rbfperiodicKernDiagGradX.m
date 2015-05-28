function gX = rbfperiodicKernDiagGradX(kern, X)

% RBFPERIODICKERNDIAGGRADX Gradient of RBFPERIODIC kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = RBFPERIODICKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the RBF derived periodic kernel matrix with respect to
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
%	RBFPERIODICKERNPARAMINIT, KERNDIAGGRADX, RBFPERIODICKERNGRADX


%	Copyright (c) 2007 Neil D. Lawrence



gX = zeros(size(X));
