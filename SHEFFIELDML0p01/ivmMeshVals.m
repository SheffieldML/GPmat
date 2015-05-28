function [X, Y, Z, varZ] = ivmMeshVals(model, limx, limy, number)

% IVMMESHVALS Give the output of the IVM for contour plot display.
%
%	Description:
%
%	[X, Y, Z] = IVMMESHVALS(MODEL, LIMX, LIMY, NUMBER) gives the output
%	of the IVM for plotting as a contour plot.
%	 Returns:
%	  X - the matrix of X values for the contour plot.
%	  Y - the matrix of Y values for the contour plot.
%	  Z - the matrix of Z values for the contour plot.
%	 Arguments:
%	  MODEL - the model to be plotted.
%	  LIMX - the x limits of the contour plot.
%	  LIMY - the y limits of the contour plot.
%	  NUMBER - the number of points to use along each side of the grid
%	   for the contour plot.
%	
%
%	See also
%	IVMCONTOUR, NCNMCONTOUR


%	Copyright (c) 2005, 2004 Neil D. Lawrence


x = linspace(limx(1), limx(2), number);
y = linspace(limy(1), limy(2), number);
[X, Y] = meshgrid(x, y);

[Z, varZ] = ivmPosteriorMeanVar(model, [X(:) Y(:)]);
Z = reshape(Z, size(X));
varZ = reshape(varZ, size(X));