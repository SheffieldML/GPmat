function K = velotransKernCompute(kern, varargin)

% VELOTRANSKERNCOMPUTE Compute the VELOTRANS kernel given the parameters and X.
%
%	Description:
%
%	K = VELOTRANSKERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the velocity translate kernel given inputs associated with rows
%	and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel. The
%	   last column of the input matrix is time.
%	  X2 - the input matrix associated with the columns of the kernel.
%	   The last column of the input matrix is time.
%
%	K = VELOTRANSKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	velocity translate kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix. The last
%	   column of the input is time.
%	
%
%	See also
%	VELOTRANSKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, VELOTRANSKERNDIAGCOMPUTE, TRANSLATEKERNCOMPUTE


%	Copyright (c) 2011 Neil D. Lawrence


  
for i = 1:length(varargin)
  t = varargin{i}(:, end);
  varargin{i}(:, end) = [];
  varargin{i} = varargin{i} - t*kern.velocity;
end
K = cmpndKernCompute(kern, varargin{:});