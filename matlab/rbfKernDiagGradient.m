function g = rbfKernDiagGradient(kern, x, covDiag)

% RBFKERNDIAGGRADIENT Compute the gradient of the RBF kernel's diagonal wrt to parameters.

% KERN

g = zeros(1, kern.nParams);
g(1) = 0;
g(2) = sum(covDiag);
