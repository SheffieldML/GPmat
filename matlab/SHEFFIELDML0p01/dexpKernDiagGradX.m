function gX = dexpKernDiagGradX(kern, x)

% DEXPKERNDIAGGRADX Gradient of the double exponential kernel's diagonal
%
%	Description:
%	with respect to x.
%	
%
%	GX = DEXPKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the double exponential kernel matrix with respect to the
%	elements of the inputs in x.
%	 Returns:
%	  GX - the gradients of the diagonal with respect to each element of
%	   x. The returned vector has the same dimensions as x.
%	 Arguments:
%	  KERN - the kernel structure for which gradients are being
%	   computed.
%	  X - the input data in the form of a design matrix.
%	
%
%	See also
%	DEXPKERNPARAMINIT, KERNDIAGGRADX, DEXPKERNGRADX


%	Copyright (c) 2009 David Luengo



gX = zeros(size(x));
