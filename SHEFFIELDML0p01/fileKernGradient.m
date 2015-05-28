function g = fileKernGradient(kern, x, varargin)

% FILEKERNGRADIENT Gradient of FILE kernel's parameters.
%
%	Description:
%
%	G = FILEKERNGRADIENT(KERN, X, PARTIAL) computes the gradient of
%	functions with respect to the stored file kernel's parameters. As
%	well as the kernel structure and the input positions, the user
%	provides a matrix PARTIAL which gives the partial derivatives of the
%	function with respect to the relevant elements of the kernel matrix.
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
%	G = FILEKERNGRADIENT(KERN, X1, X2, PARTIAL) computes the derivatives
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
%	% SEEALSO FILEKERNPARAMINIT, KERNGRADIENT, FILEKERNDIAGGRADIENT, KERNGRADX


%	Copyright (c) 2005, 2006 Neil D. Lawrence



% The last argument is covGrad
if nargin < 4
  [k, sk] = fileKernCompute(kern, x);
else
  [k, sk] = fileKernCompute(kern, x, varargin{1});
end
g(1) =  sum(sum(varargin{end}.*sk));
