function gX = sqexpKernDiagGradX(kern, x)

% SQEXPKERNDIAGGRADX Gradient of squared exponential kernel's diagonal with respect to a point x.

% IVM

gX = zeros(size(x));
