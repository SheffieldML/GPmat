function gX = lmcKernDiagGradX(kern, X)

% LMCKERNDIAGGRADX Gradient of LMC kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = LMCKERNDIAGGRADX(KERN, X) computes the gradient of the diagonal
%	of the LMC kernel matrix with respect to the elements of the design
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
%	See also
%	LMCKERNPARAMINIT, KERNDIAGGRADX


%	Copyright (c) 2010 Mauricio A. Alvarez

  
gX = zeros([kern.nout*size(X, 1) size(X,2)]);

