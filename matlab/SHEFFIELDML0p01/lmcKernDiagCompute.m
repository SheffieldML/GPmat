function k = lmcKernDiagCompute(kern, X)

% LMCKERNDIAGCOMPUTE Compute the diagonal of the LMC kernel.
%
%	Description:
%
%	K = LMCKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the kernel
%	matrix for the LMC kernel function given X.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	LMCKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, LMCKERNCOMPUTE


%	Copyright (c) 2010 Mauricio A. Alvarez


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

