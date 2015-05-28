function K = translateKernCompute(kern, varargin)

% TRANSLATEKERNCOMPUTE Compute the TRANSLATE kernel given the parameters and X.
%
%	Description:
%
%	K = TRANSLATEKERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the input space translation kernel given inputs associated with
%	rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = TRANSLATEKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	input space translation kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	cmpndKernCompute, translateKernDiagCompute
%	
%
%	See also
%	TRANSLATEKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, 


%	Copyright (c) 2007 Neil D. Lawrence


for i = 1:length(varargin)
  varargin{i} = varargin{i} - repmat(kern.centre, size(varargin{i}, 1), 1);
end
K = cmpndKernCompute(kern, varargin{:});