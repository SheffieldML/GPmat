function g = lfmKernGradient(kern, t, varargin)

% LFMKERNGRADIENT Gradient of LFM kernel's parameters.
%
%	Description:
%
%	G = LFMKERNGRADIENT(KERN, T, PARTIAL) computes the gradient of
%	functions with respect to the latent force model kernel's
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
%	   the number of rows in t.
%
%	G = LFMKERNGRADIENT(KERN, T1, T2, PARTIAL) computes the derivatives
%	as above, but input locations are now provided in two matrices
%	associated with rows and columns of the kernel matrix.
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
%	   the same number of rows as t1 and the same number of columns as t2
%	   has rows.
%	lfmXlfmKernGradient
%	
%
%	See also
%	% SEEALSO LFMKERNPARAMINIT, KERNGRADIENT, LFMKERNDIAGGRADIENT, KERNGRADX, 


%	Copyright (c) 2007 Neil D. Lawrence


if length(varargin)<2
  t2 = t;
else
  t2 = varargin{1};
end

[g1, g2] = lfmXlfmKernGradient(kern, kern, t, t2, varargin{end});

g = real(g1 + g2);


%g = real(g1);
