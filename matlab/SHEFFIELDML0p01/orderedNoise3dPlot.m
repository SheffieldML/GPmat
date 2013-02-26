function h = orderedNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% ORDEREDNOISE3DPLOT Draws a 3D or contour plot for the ORDERED noise model.
%
%	Description:
%
%	H = ORDEREDNOISE3DPLOT(NOISE, PLOTTYPE, X, Y, MU, VARSIGMA, ...)
%	draws a 3D or contour plot for the ordered categorical noise model.
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
%	ORDEREDNOISEPARAMINIT, NOISE3DPLOT, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



CZ = (CZ+noise.bias)./sqrt(CZVar);
fhandle = str2func(plotType);
fhandle(CX, CY, CZ, varargin{:});
hold on
h = [];
for i = 2:noise.C-1
  CZ = CZ - noise.widths(i-1)./sqrt(CZVar);
  h = [h; fhandle(CX, CY, CZ, varargin{:})];
end
