function h = noise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% NOISE3DPLOT Draw a 3D or contour plot for the relevant noise model.
%
%	Description:
%
%	H = NOISE3DPLOT(NOISE, PLOTTYPE, X, Y, MU, VARSIGMA, ...) draws a 3D
%	or contour plot for the relevant noise model.
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
%	NOISEPARAMINIT, NOISE3DPLOT, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence


functionName = [noise.type 'Noise3dPlot'];
if exist(functionName) == 2
  fhandle = str2func(functionName);
  h = fhandle(noise, plotType, CX, CY, CZ, CZVar, varargin{:});
end

