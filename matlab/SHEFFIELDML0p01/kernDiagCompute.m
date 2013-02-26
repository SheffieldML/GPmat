function k = kernDiagCompute(kern, x)

% KERNDIAGCOMPUTE Compute the kernel given the parameters and X.
%
%	Description:
%
%	K = KERNDIAGCOMPUTE(KERN, X) computes the diagonal of a kernel
%	matrix for the given kernel.
%	 Returns:
%	  K - vector containing computed diagonal elements of the kernel
%	   structure.
%	 Arguments:
%	  KERN - kernel structure to be computed.
%	  X - input data matrix (rows are data points) to the kernel
%	   computation.
%	
%
%	See also
%	KERNCOMPUTE


%	Copyright (c) 2003, 2004, 2005, 2006 Neil D. Lawrence


fhandle = str2func([kern.type 'KernDiagCompute']);
k = fhandle(kern, x);
