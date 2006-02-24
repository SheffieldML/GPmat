function g = rbfardKernDiagGradient(kern, x, covDiag)

% RBFARDKERNDIAGGRADIENT Compute the gradient of the RBFARD kernel's diagonal wrt to parameters.

% KERN

g = zeros(1, size(x, 2)+2);
g(2) = sum(covDiag);