function [X, Y, Z, varZ] = ivmMeshVals(model, limx, limy, number)

% IVMMESHVALS Give the output of the IVM for contour plot display.
% FORMAT
% DESC gives the output of the IVM for plotting as a contour plot.
% ARG model : the model to be plotted.
% ARG limx : the x limits of the contour plot.
% ARG limy : the y limits of the contour plot.
% ARG number : the number of points to use along each side of the
% grid for the contour plot.
% RETURN X : the matrix of X values for the contour plot.
% RETURN Y : the matrix of Y values for the contour plot.
% RETURN Z : the matrix of Z values for the contour plot.
%
% SEEALSO : ivmContour, ncnmContour
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2004

% IVM

x = linspace(limx(1), limx(2), number);
y = linspace(limy(1), limy(2), number);
[X, Y] = meshgrid(x, y);

[Z, varZ] = ivmPosteriorMeanVar(model, [X(:) Y(:)]);
Z = reshape(Z, size(X));
varZ = reshape(varZ, size(X));
