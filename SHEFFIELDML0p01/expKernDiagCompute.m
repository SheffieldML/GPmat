function k = expKernDiagCompute(kern, x)

% EXPKERNDIAGCOMPUTE Compute diagonal of EXP kernel.
%
%	Description:
%
%	K = EXPKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the kernel
%	matrix for the exponentiated kernel given a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	EXPKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, EXPKERNCOMPUTE


%	Copyright (c) 2006 Neil D. Lawrence


if kern.isStationary
  diagVal = kernDiagCompute(kern.argument, x(1, :));
  k = repmat(kern.variance*exp(diagVal)*(exp(diagVal) - 1), size(x, 1), 1);
else
  diagVal = kernDiagCompute(kern.argument, x);
  k = kern.variance*exp(diagVals).*(exp(diagVals)-1);
end
