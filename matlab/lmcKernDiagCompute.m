function k = lmcKernDiagCompute(kern, X)

% LMCKERNDIAGCOMPUTE Compute the diagonal of the LMC kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the LMC kernel 
% function given X.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : input data matrix in the form of a design matrix.
%
% SEEALSO : lmcKernParamInit, kernCompute, kernCreate, lmcKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez,  2010

% KERN

fhandle = str2func([kern.basicKernelType 'KernDiagCompute']);
if iscell(X)
    if length(X) ~= kern.nout
        error('Time information is not matched among blocks!');
    end
    dim1 = zeros(1,kern.nout);
    for i = 1:kern.nout
        dim1(i) = size(X{i}, 1);
    end
    k = zeros(sum(dim1),1);
    startOne = 1;
    endOne = 0;
    for i=1:kern.nout
        endOne = endOne + dim1(i);
        k(startOne:endOne) = kern.B(i,i)*fhandle(kern, X{i});
        startOne = endOne + 1;
    end
else
    basick = fhandle(kern, X);
    k = zeros(kern.nout*size(X,1),1);
    startOne = 1;
    endOne = 0;
    for i=1:kern.nout
        endOne = endOne + size(X,1);
        k(startOne:endOne) = kern.B(i,i)*basick;
        startOne = endOne + 1;
    end
end

