function g = whiteKernDiagGradient(kern, x, covDiag)

% WHITEKERNDIAGGRADIENT Compute the gradient of the white kernel's diagonal wrt to parameters.

% KERN

g(1) = sum(covDiag);
