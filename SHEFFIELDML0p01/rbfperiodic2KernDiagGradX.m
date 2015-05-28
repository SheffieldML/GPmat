function gX = rbfperiodic2KernDiagGradX(kern, X)

% RBFPERIODIC2KERNDIAGGRADX Gradient of RBFPERIODIC2 kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = RBFPERIODIC2KERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the RBF periodic covariance with variying period kernel
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
%
%	See also
%	RBFPERIODIC2KERNPARAMINIT, KERNDIAGGRADX, RBFPERIODIC2KERNGRADX


%	Copyright (c) 2007, 2009 Neil D. Lawrence


%	With modifications by Andreas C. Damianou 2011


%	With modifications by Michalis K. Titsias 2011



gX = zeros(size(X));
