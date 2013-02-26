function [k, sk, n2] = rbfperiodic2KernCompute(kern, x, x2)

% RBFPERIODIC2KERNCOMPUTE Compute the RBFPERIODIC2 kernel given the parameters and X.
%
%	Description:
%
%	K = RBFPERIODIC2KERNCOMPUTE(KERN, X, X2) computes the kernel
%	parameters for the RBF periodic covariance with variying period
%	kernel given inputs associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = RBFPERIODIC2KERNCOMPUTE(KERN, X) computes the kernel matrix for
%	the RBF periodic covariance with variying period kernel given a
%	design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%	
%
%	See also
%	RBFPERIODIC2KERNPARAMINIT, KERNCOMPUTE, KERNCREATE, RBFPERIODIC2KERNDIAGCOMPUTE


%	Copyright (c) 2007, 2009 Neil D. Lawrence


%	With modifications by Andreas C. Damianou 2011


%	With modifications by Michalis K. Titsias 2011


factor = kern.factor; % Default (if period is fixed: 2*pi/kern.period)
if nargin < 3
  n2 = sin(0.5*factor*(repmat(x, 1, size(x, 1)) - repmat(x', size(x, 1), 1)));
  n2 = n2.*n2;
  wi2 = (2 .* kern.inverseWidth);
  sk = exp(-n2*wi2);
else
  n2 = sin(0.5*factor*(repmat(x, 1, size(x2, 1)) - repmat(x2', size(x, 1), 1)));  
  n2 = n2.*n2;
  wi2 = (2 .* kern.inverseWidth);
  sk = exp(-n2*wi2);
end
k = kern.variance*sk;
% Test kernel with: kernTest('rbfperiodic2',1)

  