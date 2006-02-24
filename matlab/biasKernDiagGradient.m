function g = biasKernDiagGradient(kern, x, covDiag)

% BIASKERNDIAGGRADIENT Compute the gradient of the bias kernel's diagonal wrt to parameters.

% KERN

g(1) = sum(covDiag);
