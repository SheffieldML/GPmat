function gX = linKernGradX(kern, x, X2)

% LINKERNGRADX Gradient of linear kernel with respect to a point X.

% KERN

% KERN


gX = kern.variance.*X2;