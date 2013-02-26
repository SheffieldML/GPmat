function K = whiteblockKernCompute(kern, x, x2)

% WHITEBLOCKKERNCOMPUTE Compute the WHITEBLOCK kernel.
%
%	Description:
%
%	K = WHITEBLOCKKERNCOMPUTE(KERN, X, X2) computes the kernel matrix
%	from the white noise block kernel given inputs associated with rows
%	and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the inpute matrix associated with the columns of the kernel.
%
%	K = WHITEBLOCKKERNCOMPUTE(KERN, X) computes the kernel matrix for
%	the white noise block kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	WHITEBLOCKKERNPARAMINIT


%	Copyright (c) 2010 Mauricio A. Alvarez


if nargin < 3
  diagk = whiteblockKernDiagCompute(kern, x);  
  K = sparseDiag(diagk);
else
    if iscell(x)
        dim1 = zeros(kern.nout,1);
        dim2 = zeros(kern.nout,1);
        for i=1:kern.nout
            dim1(i) = size(x{i},1);
            dim2(i) = size(x2{i},1);
        end
        K = spalloc(sum(dim1), sum(dim2), 0);
    else
        K = spalloc(kern.nout*size(x, 1), kern.nout*size(x2, 1), 0);
    end
end
