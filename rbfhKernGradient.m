function g = rbfhKernGradient(kern, x1, x2, covGrad)

% RBFHKERNGRADIENT Gradient of RBFH kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% radial basis function heat kernel's parameters. As well as the kernel
% structure and the input positions, the user provides a matrix PARTIAL
% which gives the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations for which the gradients are being
% computed.
% ARG covGrad : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations associated with the rows of the
% kernel matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix.
% ARG covGrad : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO: rbfhKernParamInit, kernGradient
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 4
    covGrad = x2;
    x2 = x1;
end
if size(x1, 2) ~= 2 || size(x2, 2) ~= 2
    error('Input can only have two columns');
end

% Split the domain into time domain and spatial domain and account for
% missing values. If there are no missing values the computation of the
% kernel is a pointwise prodruct, otherwise it is a kronecker product.
t1 = x1(x1(:,1)~=Inf,1);
t2 = x2(x2(:,1)~=Inf,1);
s1 = x1(x1(:,2)~=Inf,2);
s2 = x2(x2(:,2)~=Inf,2);
if (length(t1) == length(s1)) && (length(t2) == length(s2))
    ut1 = unique(t1);
    ut2 = unique(t2);
    us1 = unique(s1);
    us2 = unique(s2);
    if (length(ut1)*length(us1) == length(t1)) && ...
            (length(ut2)*length(us2) == length(t2))
        t1 = ut1; s1 = us1; t2 = ut2; s2 = us2;
        isPointwise = false;
    else
        isPointwise = true;
    end
else
    isPointwise = false;
end

kern.rbf.inverseWidth = kern.inverseWidthTime;
Kt = rbfKernCompute(kern.rbf, t1, t2);
kern.rbf.inverseWidth = kern.inverseWidthSpace;
Ks = rbfKernCompute(kern.rbf, s1, s2);

if isPointwise
    covGradt = covGrad.*Ks;
    covGrads = covGrad.*Kt;
else
    covGradt = zeros(length(t1), length(t2));
    covGrads = zeros(length(s1), length(s2));
    endOne = 0;
    startOne = 1;
    for k=1:length(t1)
        endOne = endOne + length(s1);
        startTwo = 1;
        endTwo = 0;
        for l=1:length(t2)
            endTwo = endTwo + length(s2);
            covGradt(k,l) = sum(sum(covGrad(startOne:endOne, startTwo:endTwo).*Ks));
            startTwo = endTwo + 1;
        end
        startOne = endOne + 1;
    end
    for k=1:length(s1)
        indRows = (k:length(s1):(length(t1)*length(s1)))';
        for l=1:length(s2)
            indCols = l:length(s2):(length(t2)*length(s2));
            covGrads(k,l) = sum(sum(covGrad(indRows, indCols).*Kt));
        end
    end
end

kern.rbf.inverseWidth = kern.inverseWidthTime;
gt = rbfKernGradient(kern.rbf, t1, t2, covGradt);
kern.rbf.inverseWidth = kern.inverseWidthSpace;
gs = rbfKernGradient(kern.rbf, s1, s2, covGrads);

g = [gt(1) gs(1)];

