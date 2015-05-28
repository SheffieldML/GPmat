function k = invcmpndKernDiagCompute(kern, x)

% INVCMPNDKERNDIAGCOMPUTE Compute diagonal of INVCMPND kernel.
%
%	Description:
%
%	K = INVCMPNDKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the inv. precision compound kernel given a design
%	matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	INVCMPNDKERNPARAMINIT, CMPNDKERNDIAGCOMPUTE, KERNCREATE, INVCMPNDKERNCOMPUTE


%	Copyright (c) 2012 Andreas C. Damianou



% Unlike the cmpnd  version, here the nonlinearity of the inverses means
% that we cannot trivially compute the diagonal. The naive way here is to
% just compute the whole kernel and return the diagonal.

k = diag(invcmpndKernCompute(kern, x));