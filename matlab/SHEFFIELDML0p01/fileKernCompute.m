function [k, sk] = fileKernCompute(kern, varargin)

% FILEKERNCOMPUTE Compute the FILE kernel given the parameters and X.
%
%	Description:
%
%	K = FILEKERNCOMPUTE(KERN, INDEX1, INDEX2) computes the kernel
%	parameters for the stored file kernel given inputs associated with
%	rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  INDEX1 - the row indices of the kernel matrix to return.
%	  INDEX2 - the column indices of the kernel matrix to return.
%
%	K = FILEKERNCOMPUTE(KERN, INDEX) computes the kernel matrix for the
%	stored file kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  INDEX - indices of the kernel matrix to return.
%	
%
%	See also
%	FILEKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, FILEKERNDIAGCOMPUTE


%	Copyright (c) 2005, 2006 Neil D. Lawrence


sk = fileKernRead(kern, varargin{:});
k = kern.variance*sk;
