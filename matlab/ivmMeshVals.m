function [X, Y, Z, varZ] = ivmMeshVals(model, limx, limy, number)

% IVMMESHVALS Give the output of the IVM for contour plot display.

% IVM

x = linspace(limx(1), limx(2), number);
y = linspace(limy(1), limy(2), number);
[X, Y] = meshgrid(x, y);

[Z, varZ] = ivmPosteriorMeanVar(model, [X(:) Y(:)]);
Z = reshape(Z, size(X));
varZ = reshape(varZ, size(X));