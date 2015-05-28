function gX = lfmKernDiagGradX(kern, X)

% LFMKERNDIAGGRADX Gradient of LFM kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = LFMKERNDIAGGRADX(KERN, X) computes the gradient of the diagonal
%	of the single input motif kernel matrix with respect to the elements
%	of the design matrix given in X.
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
%	LFMKERNPARAMINIT, KERNDIAGGRADX, LFMKERNGRADX


%	Copyright (c) 2006 Neil D. Lawrence


error('lfmKernDiagGradX not yet implemented.');