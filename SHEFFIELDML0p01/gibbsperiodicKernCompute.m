function [K, sk, n2, w2, l, l2] = gibbsperiodicKernCompute(kern, x, x2)

% GIBBSPERIODICKERNCOMPUTE Compute the GIBBSPERIODIC kernel given the parameters and X.
%
%	Description:
%
%	K = GIBBSPERIODICKERNCOMPUTE(KERN, X, X2) computes the kernel
%	parameters for the Gibbs-kernel derived periodic kernel given inputs
%	associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = GIBBSPERIODICKERNCOMPUTE(KERN, X) computes the kernel matrix for
%	the Gibbs-kernel derived periodic kernel given a design matrix of
%	inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	GIBBSPERIODICKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, GIBBSPERIODICKERNDIAGCOMPUTE


%	Copyright (c) 2007, 2009 Neil D. Lawrence



fhandle = str2func([kern.lengthScaleTransform, 'Transform']);
l = fhandle(modelOut(kern.lengthScaleFunc, x), 'atox');
if nargin < 3
  n2 = sin(0.5*(repmat(x, 1, size(x, 1)) - repmat(x', size(x, 1), 1)));  
  n2 = 4*n2.*n2;
  L = repmat(l.*l, 1, size(l, 1));
  w2 = L + L';
  sk = ((2*l*l')./w2).*exp(-n2./w2);
else
  n2 = sin(0.5*(repmat(x, 1, size(x2, 1)) - repmat(x2', size(x, 1), 1)));  
  n2 = 4*n2.*n2;
  l2 = fhandle(modelOut(kern.lengthScaleFunc, x2), 'atox');
  w2 = repmat(l.*l, 1, size(l2, 1))+repmat(l2.*l2, 1, size(l, 1))';
  sk = ((2*l*l2')./w2).*exp(-n2./w2);
end
K = kern.variance*sk;
