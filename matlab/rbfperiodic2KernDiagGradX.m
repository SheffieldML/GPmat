function gX = rbfperiodic2KernDiagGradX(kern, X)

% RBFPERIODIC2KERNDIAGGRADX Gradient of RBFPERIODIC2 kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = RBFPERIODIC2KERNDIAGGRADX(KERN, X) computes the gradient of the
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
%	RBFPERIODIC2KERNPARAMINIT, KERNDIAGGRADX, RBFPERIODIC2KERNGRADX


%	Copyright (c) 2007 Neil D. Lawrence
% 	rbfperiodicKernDiagGradX.m CVS version 1.1
% 	rbfperiodicKernDiagGradX.m SVN version 1
% 	last update 2007-02-03T09:25:10.000000Z


gX = zeros(size(X));
