function gX = linardKernGradX(kern, x, X2)

% LINARDKERNGRADX Gradient of linear ARD kernel with respect to a point x.

% KERN


scales = sparse(diag(kern.inputScales));
X2 = X2*scales;

gX = kern.variance.*X2;