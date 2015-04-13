function KPart = kernelDiffX(X, theta, q);

% KERNELDIFFX Compute the gradient of the kernel with respect to X.

% Since the result is all zeros apart from one row/column combination,
% and the result is symetric, here we return a matrix whose rows are the
% vectors which give the row/column combinations  

theta = thetaConstrain(theta);

latentDim = size(X, 2);
numData = size(X, 1);

K1 = repmat(X(:, q), 1, numData);
K2 = repmat(X(:, q)', numData, 1);
KPart =  -theta(1)*(K1-K2).*kernel(X, X, theta);
