function K = lmcKernCompute(kern, X, X2)

% LMCKERNCOMPUTE Compute the LMC kernel given the parameters and X.
% FORMAT
% DESC computes the kernel matrix for the linear model of coregionalization
% kernel given inputs associated with X (rows) and X2 (columns).
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 the input matrix associated with the columns of the kernel.
%
% FORMAT
% DESC computes the kernel matrix for the LMC kernel function given X.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : input data matrix in the form of a design matrix.
%
% SEEALSO : lmcKernParamInit, kernCompute, kernCreate
%
% COPYRIGHT : Mauricio A. Alvarez,  2010

% KERN


fhandle = str2func([kern.basicKernelType 'KernCompute']);

if iscell(X)
    if length(X) ~= kern.nout
        error('Time information is not matched among blocks!');
    end
    dim1 = zeros(1,kern.nout);
    dim2 = zeros(1,kern.nout);
    for i = 1:kern.nout
        dim1(i) = size(X{i}, 1);
        if nargin>2
            if length(X) ~= length(X2)
                error('Time information is not matched within the block!');
            end
            dim2(i) = size(X2{i}, 1);
        else
            dim2(i) = dim1(i);
        end
    end
    K = zeros(sum(dim1), sum(dim2));
    for i = 1:kern.nout
        startOne = sum(dim1(1:(i-1))) + 1;
        endOne = sum(dim1(1:i));
        startThree = sum(dim2(1:(i-1))) + 1;
        endThree = sum(dim2(1:i));
        if nargin<3
            K(startOne:endOne, startThree:endThree) = kern.B(i,i)*fhandle(kern, X{i});
        else
            K(startOne:endOne, startThree:endThree) = kern.B(i,i)*fhandle(kern, X{i}, X2{i});
        end
        for j = 1:i-1
            startTwo = sum(dim2(1:(j-1))) + 1;
            endTwo =  sum(dim2(1:j));

            if nargin<3
                K(startOne:endOne, startTwo:endTwo) = kern.B(i,j)*fhandle(kern, X{i}, X{j});
                K(startTwo:endTwo, startOne:endOne) = K(startOne:endOne, ...
                    startTwo:endTwo)';
            else
                K(startOne:endOne, startTwo:endTwo) = kern.B(i,j)*fhandle(kern, X{i}, X2{j});                
                startFour = sum(dim1(1:(j-1))) + 1;
                endFour =  sum(dim1(1:j));
                K(startFour:endFour, startThree:endThree) = (kern.B(j,i)*fhandle(kern, X2{i}, X{j})');
            end
        end
    end
else
    if nargin < 3
        X2 = X;
    end
    basicK = fhandle(kern, X, X2);
    K = zeros(kern.nout*size(X,1), kern.nout*size(X2,1));

    % This is a basic implementation of the Kronecker product. Most
    % computational efficient implementations can be found in signal
    % processing literature.
    startOne = 1;
    endOne = 0;
    startThree = 1;
    endThree = 0;
    for i=1:kern.nout
        endOne = endOne + size(X,1);
        endThree = endThree + size(X2,1);
        K(startOne:endOne, startThree:endThree) = kern.B(i,i)*basicK;
        startTwo = 1;
        endTwo = 0;
        startFour = 1;
        endFour = 0;
        for j=1:i-1
            endTwo = endTwo + size(X2,1);
            K(startOne:endOne, startTwo:endTwo) = kern.B(i,j)*basicK;
            if nargin < 3
                K(startTwo:endTwo, startOne:endOne) = kern.B(j,i)*(basicK');
            else
                endFour = endFour + size(X,1);
                K(startFour:endFour, startThree:endThree) = kern.B(j,i)*basicK;
                startFour = endFour + 1;
            end
            startTwo = endTwo + 1;
        end
        startOne = endOne + 1;
        startThree = endThree + 1;
    end
end
