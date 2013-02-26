function gX = ratquadKernDiagGradX(kern, X)

% RATQUADKERNDIAGGRADX Gradient of RATQUAD kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = RATQUADKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the rational quadratic kernel matrix with respect to the
%	elements of the design matrix given in X.
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
%	RATQUADKERNPARAMINIT, KERNDIAGGRADX, RATQUADKERNGRADX


%	Copyright (c) 2006 Neil D. Lawrence


gX = zeros(size(X));
