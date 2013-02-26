function gX = rbfhKernDiagGradX(kern, X)

% RBFHKERNDIAGGRADX Gradient of RBFH kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = RBFHKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the radial basis function heat kernel matrix with
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
%	RBFHKERNPARAMINIT, KERNDIAGGRADX


%	Copyright (c) 2010 Mauricio A. Alvarez


gX = zeros(size(X));

