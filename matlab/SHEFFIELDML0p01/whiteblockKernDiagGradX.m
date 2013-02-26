function gX = whiteblockKernDiagGradX(kern, X)

% WHITEBLOCKKERNDIAGGRADX Gradient of WHITEBLOCK kernel's diagonal wrt X.
%
%	Description:
%
%	GX = WHITEBLOCKKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the white noise block kernel matrix with respect to the
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
%	WHITEBLOCKKERNPARAMINIT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



gX = zeros([kern.nout*size(X,1) size(X,2)]);
