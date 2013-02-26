function [K, sk, n2, w2, l, l2] = gibbsKernCompute(kern, x, x2)

% GIBBSKERNCOMPUTE Compute the GIBBS kernel given the parameters and X.
%
%	Description:
%
%	K = GIBBSKERNCOMPUTE(KERN, X, X2) computes the kernel parameters for
%	Mark Gibbs's non-stationary kernel given inputs associated with rows
%	and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = GIBBSKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	Mark Gibbs's non-stationary kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	GIBBSKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, GIBBSKERNDIAGCOMPUTE


%	Copyright (c) 2006, 2009 Neil D. Lawrence



fhandle = str2func([kern.lengthScaleTransform, 'Transform']);
l = fhandle(modelOut(kern.lengthScaleFunc, x), 'atox');
if nargin < 3
  n2 = dist2(x, x);
  L = repmat(l.*l, 1, size(l, 1));
  w2 = L + L';
  sk = ((2*l*l')./w2).^(kern.inputDimension/2).*exp(-n2./w2);
  K = kern.variance*sk;
else
  n2 = dist2(x, x2);
  l2 = fhandle(modelOut(kern.lengthScaleFunc, x2), 'atox');
  w2 = repmat(l.*l, 1, size(l2, 1))+repmat(l2.*l2, 1, size(l, 1))';
  sk = ((2*l*l2')./w2).^(kern.inputDimension/2).*exp(-n2./w2);
  K = kern.variance*sk;
end
