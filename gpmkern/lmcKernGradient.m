function g = lmcKernGradient(kern, X, X2, covGrad)

% LMCKERNGRADIENT Gradient of LMC kernel's parameters.
% FORMAT
% DESC computes the gradient of parameters associated to the LMC kernel.
% As well as the kernel structure and the input positions, the user
% provides a matrix COVGRAD which gives the partial derivatives of the
% function with respect to the relevant elements of the kernel matrix.
% RETURN g:  gradients of the function of interest with respect to the
% kernel parameters. The ordering of the vector should match that
% provided by the function kernExtractParam.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG X : the input locations for which the gradients are being computed.
% ARG covGrad : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes the
% form of a square matrix of dimension  numData, where numData is the
% number of rows in X.
%
% FORMAT
% DESC  computes the derivatives as above, but input locations are now
% provided in two matrices associated with rows and columns of the kernel
% matrix.
% RETURN g : gradients of the function of interest with respect to the
% kernel parameters.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG X : the input locations associated with the rows of the kernel matrix.
% ARG X2 : the input locations associated with the columns of the kernel
% matrix.
% ARG covGrad : matrix of partial derivatives of the function of interest
% with respect to the kernel matrix. The matrix should have the same number
% of rows as X1 and the same number of columns as X2 has rows.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010.

% KERN

fhandle = str2func([kern.basicKernelType 'KernCompute']);
fhandleGrad = str2func([kern.basicKernelType 'KernGradient']);
gBK = zeros(1, kern.nParamsBK);
gPartialB = zeros(kern.nout);

if iscell(X)
    if nargin>3 && ~iscell(X2)
        error('Time course information is not matched in Cell format!');
    end
    % Collate arguments.
    dim1 = zeros(1, kern.nout);
    dim2 = zeros(1, kern.nout);
    for i=1:kern.nout
        dim1(i) = size(X{i}, 1);
        if nargin > 3
            dim2(i) = size(X2{i}, 1);
        else
            dim2(i) = dim1(i);
            covGrad = X2;
        end
    end
    for i = 1:kern.nout
        startOne = sum(dim1(1:(i-1)))+1;
        endOne = sum(dim1(1:i));
        startThree = sum(dim2(1:(i-1))) + 1;
        endThree = sum(dim2(1:i));
        if nargin > 3
            gBK = gBK + kern.B(i,i)*fhandleGrad(kern, X{i}, X2{i}, covGrad(startOne:endOne, ...
                startThree:endThree));
            basicK = fhandle(kern, X{i}, X2{i});
        else
            gBK = gBK + kern.B(i,i)*fhandleGrad(kern, X{i}, covGrad(startOne:endOne,...
                startThree:endThree));
            basicK = fhandle(kern, X{i});
        end
        gPartialB(i,i) = sum(sum(covGrad(startOne:endOne, startThree:endThree).*basicK));
        for j = 1:i-1
            startTwo = sum(dim2(1:(j-1))) + 1;
            endTwo =  sum(dim2(1:j));
            if nargin > 3
                g2  = kern.B(i,j)*fhandleGrad(kern, X{i}, X2{j}, covGrad(startOne:endOne, ...
                    startTwo:endTwo));
                basicK = fhandle(kern, X{i}, X2{j});
            else
                g2  = kern.B(i,j)*fhandleGrad(kern, X{i}, X{j}, covGrad(startOne:endOne, ...
                    startTwo:endTwo));
                basicK = fhandle(kern, X{i}, X{j});
            end
            gBK = gBK + 2*g2;
            gPartialB(i,j) = sum(sum(covGrad(startOne:endOne, startTwo:endTwo).*basicK));
            gPartialB(j,i) = gPartialB(i,j);
        end
    end
else
    if nargin < 4
        covGrad = X2;
        X2 = X;
    end
    basicK = fhandle(kern, X, X2);
    startOne = 1;
    endOne = 0;
    startThree = 1;
    endThree = 0;
    for i=1:kern.nout
        endOne = endOne + size(X,1);
        endThree = endThree + size(X2,1);
        gBK = gBK + kern.B(i,i)*fhandleGrad(kern, X, X2, ...
            covGrad(startOne:endOne, startThree:endThree));
        gPartialB(i,i) = sum(sum(covGrad(startOne:endOne, startThree:endThree).*basicK));
        startTwo = 1;
        endTwo = 0;
        startFour = 1;
        endFour = 0;
        for j=1:i-1
            endTwo = endTwo + size(X2,1);
            g2  = kern.B(i,j)*fhandleGrad(kern, X, X2, covGrad(startOne:endOne, ...
                startTwo:endTwo));
            gBK = gBK + g2;
            gPartialB(i,j) = sum(sum(covGrad(startOne:endOne, startTwo:endTwo).*basicK));
            if nargin < 3
                gBK = gBK + g2;
                gPartialB(j,i) = gPartialB(i,j);
            else
                endFour = endFour + size(X,1);
                g3 = kern.B(j,i)*fhandleGrad(kern, X, X2, covGrad(startFour:endFour, ...
                    startThree:endThree));
                gBK = gBK + g3;
                gPartialB(j,i) = sum(sum(covGrad(startFour:endFour, startThree:endThree).*basicK));
                startFour = endFour + 1;
            end
            startTwo = endTwo + 1;
        end
        startOne = endOne + 1;
        startThree = endThree + 1;
    end
end

Jij = zeros(kern.nout, kern.rankCorregMatrix);
Jji = zeros(kern.rankCorregMatrix, kern.nout);
gB = zeros(kern.nout, kern.rankCorregMatrix);
for i=1:kern.nout
    for j=1:kern.rankCorregMatrix
        Jij(i,j) = 1; Jji(j,i) = 1;
        gB(i,j) = sum(sum((kern.A*Jij' + Jji'*kern.A').*gPartialB));
        Jij(i,j) = 0; Jji(j,i) = 0;
    end
end

g = [gBK gB(:)'];

