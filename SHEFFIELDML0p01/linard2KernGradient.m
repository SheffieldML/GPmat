function g = linard2KernGradient(kern, x, varargin)

% LINARD2KERNGRADIENT Gradient of LINARD2 kernel's parameters.
%
%	Description:
%
%	G = LINARD2KERNGRADIENT(KERN, X, PARTIAL) computes the gradient of
%	functions with respect to the automatic relevance determination
%	linear kernel's parameters. As well as the kernel structure and the
%	input positions, the user provides a matrix PARTIAL which gives the
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
%	G = LINARD2KERNGRADIENT(KERN, X1, X2, PARTIAL) computes the
%	derivatives as above, but input locations are now provided in two
%	matrices associated with rows and columns of the kernel matrix.
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
%
%	See also
%	% SEEALSO LINARD2KERNPARAMINIT, KERNGRADIENT, LINARD2KERNDIAGGRADIENT, KERNGRADX


%	Copyright (c) 2004, 2005, 2006, 2009 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias



  g = zeros(1, size(x, 2));
  if nargin < 4
    [k, sk] = linard2KernCompute(kern, x);
  else
    [k, sk] = linard2KernCompute(kern, x, varargin{1});
  end
  if nargin < 4
    for i = 1:size(x, 2)
      g(i) =  x(:, i)'*varargin{end}*x(:, i);
    end
  else
    for i = 1:size(x, 2)
      g(i) =  x(:, i)'*varargin{end}*varargin{1}(:, i);
    end
  end
end