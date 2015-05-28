function gX = diagKernDiagGradX(kern, X)

% DIAGKERNDIAGGRADX Gradient of DIAG kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = DIAGKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the diagonal noise covariance function kernel matrix
%	with respect to the elements of the design matrix given in X.
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
%	DIAGKERNPARAMINIT, KERNDIAGGRADX, DIAGKERNGRADX


%	Copyright (c) 2011 Neil D. Lawrence


  if size(X, 2)>1
    error('Diag kernel requires 1-dimensional input.')
  end
  trans = str2func([kern.trans, 'Transform']);
  vars = trans(X, 'atox');
  gX = kern.variance*trans(vars, 'gradfact');
end