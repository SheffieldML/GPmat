function k = whiteblockKernDiagCompute(kern, x)

% WHITEBLOCKKERNDIAGCOMPUTE Compute diagonal of WHITEBLOCK kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the white noise block
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : whiteblockKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if iscell(x)
    dim1 = zeros(1, kern.nout);
    for i =1:kern.nout
       dim1(i) = size(x{i},1); 
    end
    k = zeros(sum(dim1),1);
    startOne = 1;
    endOne = 0;
    for i=1:kern.nout
        endOne = endOne + dim1(i);
        k(startOne:endOne) = kern.variance(i);
        startOne = endOne + 1;
    end    
else
    k = zeros(kern.nout*size(x,1),1);
    startOne = 1;
    endOne = 0;
    for i=1:kern.nout
        endOne = endOne + size(x,1);
        k(startOne:endOne) = kern.variance(i);
        startOne = endOne + 1;
    end
end    
