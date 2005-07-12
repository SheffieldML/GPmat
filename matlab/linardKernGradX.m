function gX = linardKernGradX(kern, X, X2)

% LINARDKERNGRADX Gradient of linear ARD kernel with respect to X.

% KERN

scales = sparse(diag(kern.inputScales));
X2 = X2*scales;

gX = repmat(kern.variance.*X2, [1 1 size(X, 1)]);
