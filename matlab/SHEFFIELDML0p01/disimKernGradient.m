function g = disimKernGradient(kern, t, varargin)

% DISIMKERNGRADIENT Gradient of DISIM kernel's parameters.
%
%	Description:
%
%	G = DISIMKERNGRADIENT(KERN, T, PARTIAL) computes the gradient of
%	functions with respect to the single input motif kernel's
%	parameters. As well as the kernel structure and the input positions,
%	the user provides a matrix PARTIAL which gives the partial
%	derivatives of the function with respect to the relevant elements of
%	the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters. The ordering of the vector should match that
%	   provided by the function kernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  T - the input locations for which the gradients are being
%	   computed.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The argument takes the
%	   form of a square matrix of dimension  numData, where numData is
%	   the number of rows in X.
%
%	G = DISIMKERNGRADIENT(KERN, T1, T2, PARTIAL) computes the
%	derivatives as above, but input locations are now provided in two
%	matrices associated with rows and columns of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  T1 - the input locations associated with the rows of the kernel
%	   matrix.
%	  T2 - the input locations associated with the columns of the kernel
%	   matrix.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The matrix should have
%	   the same number of rows as X1 and the same number of columns as X2
%	   has rows.
%	
%	
%
%	See also
%	% SEEALSO DISIMKERNPARAMINIT, KERNGRADIENT, DISIMKERNDIAGGRADIENT, KERNGRADX


%	Copyright (c) 2006 Neil D. Lawrence
%	Copyright (c) 2007 Antti Honkela


if length(varargin)<2
  t2 = t;
else
  t2 = varargin{1};
end

[g1, g2] = disimXdisimKernGradient(kern, kern, t, t2, varargin{end});

g = g1 + g2;
