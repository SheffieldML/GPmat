function k = rbfKernDiagCompute(kern, x)

% RBFKERNDIAGCOMPUTE Compute diagonal of RBF kernel.
%
%	Description:
%
%	K = RBFKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the kernel
%	matrix for the radial basis function kernel given a design matrix of
%	inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	RBFKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, RBFKERNCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


%	With modifications by Mauricio Alvarez 2009, David Luengo, 2009


if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
    k = repmat(kern.variance * sqrt(kern.inverseWidth/(2*pi)), size(x, 1), 1);
else
    k = repmat(kern.variance, size(x, 1), 1);
end
