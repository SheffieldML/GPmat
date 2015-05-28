function k = velotransKernDiagCompute(kern, x)

% VELOTRANSKERNDIAGCOMPUTE Compute diagonal of VELOTRANS kernel.
%
%	Description:
%
%	K = VELOTRANSKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the velocity translate kernel given a design
%	matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix. The last
%	   column of the input data matrix should be time.
%	
%
%	See also
%	VELOTRANSKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, VELOTRANSKERNCOMPUTE, TRANSLATEKERNDIAGCOMPUTE


%	Copyright (c) 2011 Neil D. Lawrence

  
t = x(:, end);
xPass = x(:, 1:end-1);
xPass = xPass - t*kern.velocity;
k = cmpndKernDiagCompute(kern, xPass);