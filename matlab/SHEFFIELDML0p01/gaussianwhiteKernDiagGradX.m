function gX = gaussianwhiteKernDiagGradX(kern, X)

% GAUSSIANWHITEKERNDIAGGRADX Gradient of gaussian white kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = GAUSSIANWHITEKERNDIAGGRADX(KERN, X) computes the gradient of
%	the diagonal of the gaussian white kernel matrix with respect to the
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
%	GAUSSIANWHITEKERNPARAMINIT, KERNDIAGGRADX, GAUSSIANWHITEKERNGRADX


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence

  
gX = zeros(size(X));