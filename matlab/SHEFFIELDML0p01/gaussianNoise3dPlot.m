function h = gaussianNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% GAUSSIANNOISE3DPLOT Draws a 3D or contour plot for the GAUSSIAN noise model.
%
%	Description:
%
%	H = GAUSSIANNOISE3DPLOT(NOISE, PLOTTYPE, X, Y, MU, VARSIGMA, ...)
%	draws a 3D or contour plot for the Gaussian noise model.
%	 Returns:
%	  H - the gradients of the diagonal with respect to each element of
%	   X. The returned matrix has the same dimensions as X.
%	 Arguments:
%	  NOISE - the noise structure for which the plot is required.
%	  PLOTTYPE - string containing the name of the plotting function
%	   (for example mesh, contour).
%	  X - the input X data in the form of a 'mesh' matrix.
%	  Y - the input Y data in the form of a 'mesh' matrix.
%	  MU - the input mean in the form of a 'mesh' matrix.
%	  VARSIGMA - the input variance in the form of a 'mesh' matrix.
%	  ... - optional additional arguments for the given plot type.
%	
%
%	See also
%	GAUSSIANNOISEPARAMINIT, NOISE3DPLOT, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



CZ = (CZ+noise.bias);
fhandle = str2func(plotType);
h = fhandle(CX, CY, CZ, varargin{:});
