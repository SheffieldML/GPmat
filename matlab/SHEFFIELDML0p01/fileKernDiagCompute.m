function k = fileKernDiagCompute(kern, x)

% FILEKERNDIAGCOMPUTE Compute diagonal of FILE kernel.
%
%	Description:
%
%	K = FILEKERNDIAGCOMPUTE(KERN, INDEX) computes the diagonal of the
%	kernel matrix for the stored file kernel given a design matrix of
%	inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  INDEX - indices of the diagonal to return.
%	
%
%	See also
%	FILEKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, FILEKERNCOMPUTE


%	Copyright (c) 2005, 2006 Neil D. Lawrence



k = kern.variance*fileKernRead(kern, x, 'diag');
