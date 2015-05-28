function K = whiteblockKernCompute(kern, x, x2)

% WHITEBLOCKKERNCOMPUTE Compute the WHITEBLOCK kernel. 
% FORMAT
% DESC computes the kernel matrix from the white noise block
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the inpute matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the white noise block
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : whiteblockKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

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
