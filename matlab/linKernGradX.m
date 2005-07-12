function gX = linKernGradX(kern, X, X2)

% LINKERNGRADX Gradient of linear kernel with respect to X.

% KERN

gX = repmat(kern.variance.*X2, [1 1 size(X, 1)]);