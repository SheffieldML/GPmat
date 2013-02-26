function g = diagKernGradient(kern, x, varargin)

% DIAGKERNGRADIENT Gradient of DIAG kernel's parameters.
%
%	Description:
%
%	G = DIAGKERNGRADIENT(KERN, X, PARTIAL) computes the gradient of
%	functions with respect to the diagonal noise covariance function
%	kernel's parameters. As well as the kernel structure and the input
%	positions, the user provides a matrix PARTIAL which gives the
%	partial derivatives of the function with respect to the relevant
%	elements of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters. The ordering of the vector should match that
%	   provided by the function kernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  X - the input locations for which the gradients are being
%	   computed.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The argument takes the
%	   form of a square matrix of dimension  numData, where numData is
%	   the number of rows in X.
%
%	G = DIAGKERNGRADIENT(KERN, X1, X2, PARTIAL) computes the derivatives
%	as above, but input locations are now provided in two matrices
%	associated with rows and columns of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  X1 - the input locations associated with the rows of the kernel
%	   matrix.
%	  X2 - the input locations associated with the columns of the kernel
%	   matrix.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The matrix should have
%	   the same number of rows as X1 and the same number of columns as X2
%	   has rows.
%	
%
%	See also
%	% SEEALSO DIAGKERNPARAMINIT, KERNGRADIENT, DIAGKERNDIAGGRADIENT, KERNGRADX


%	Copyright (c) 2011 Neil D. Lawrence


  if nargin < 4
    trans = str2func([kern.trans, 'Transform']);
    g(1, 1) = sum(diag(varargin{end}).*trans(x, 'atox'));
  else
    g = 0;
  end
end