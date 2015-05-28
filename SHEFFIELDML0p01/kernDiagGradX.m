function k = kernDiagGradX(kern, x)

% KERNDIAGGRADX Compute the gradient of the  kernel wrt X.
%
%	Description:
%
%	GX = KERNDIAGGRADX(KERN, X) computes the gradient of the diagonal of
%	the kernel matrix with respect to the elements of the design matrix
%	given in X.
%	 Returns:
%	  GX - the gradients of the diagonal with respect to each element of
%	   X. The returned matrix has the same dimensions as X.
%	 Arguments:
%	  KERN - the kernel structure for which gradients are being
%	   computed.
%	  X - the input data in the form of a design matrix.
%
%	See also
%	KERNDIAGGRADX, KERNGRADX


fhandle = str2func([kern.type 'KernDiagGradX']);
k = fhandle(kern, x);
