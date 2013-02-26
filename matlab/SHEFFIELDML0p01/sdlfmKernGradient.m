function [g, covGradLocal] = sdlfmKernGradient(kern, t, varargin)

% SDLFMKERNGRADIENT Gradient of SDLFM kernel's parameters.
%
%	Description:
%
%	[G, COVGRADLOCAL] = SDLFMKERNGRADIENT(KERN, T, PARTIAL) computes the
%	gradient of functions with respect to the switching latent force
%	model kernel's parameters. As well as the kernel structure and the
%	input positions, the user provides a matrix PARTIAL which gives the
%	partial derivatives of the function with respect to the relevant
%	elements of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters. The ordering of the vector should match that
%	   provided by the function kernExtractParam.
%	  COVGRADLOCAL - partial derivatives for the initial conditions of
%	   the first interval.
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
%	[G, COVGRADLOCAL] = SDLFMKERNGRADIENT(KERN, T1, T2, PARTIAL)
%	computes the derivatives as above, but input locations are now
%	provided in two matrices associated with rows and columns of the
%	kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters.
%	  COVGRADLOCAL - partial derivatives for the initial conditions of
%	   the first interval.
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
%	% SEEALSO SDLFMKERNPARAMINIT, KERNGRADIENT, SDLFMKERNDIAGGRADIENT, KERNGRADX, 


%	Copyright (c) 2010 Mauricio A. Alvarez


if length(varargin)<2
  t2 = t;
else
  t2 = varargin{1};
end

[g1, g2, covGradLocal] = sdlfmXsdlfmKernGradient(kern, kern, t, t2, varargin{end});

g = real(g1 + g2);
