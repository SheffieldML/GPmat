function [k, sk] = linard2KernCompute(kern, x, x2)

% LINARD2KERNCOMPUTE Compute the LINARD2 kernel given the parameters and X.
%
%	Description:
%
%	K = LINARD2KERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the automatic relevance determination linear kernel given inputs
%	associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = LINARD2KERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	automatic relevance determination linear kernel given a design
%	matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	LINARD2KERNPARAMINIT, KERNCOMPUTE, KERNCREATE, LINARD2KERNDIAGCOMPUTE


%	Copyright (c) 2004, 2005, 2006, 2009 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias



    
if nargin < 3
  sk = x*sparse(diag(kern.inputScales))*x';
else  
  sk = x*sparse(diag(kern.inputScales))*x2';
end
k = sk;