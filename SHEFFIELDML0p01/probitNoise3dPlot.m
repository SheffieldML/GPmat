function h = probitNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% PROBITNOISE3DPLOT Draws a 3D or contour plot for the PROBIT noise model.
%
%	Description:
%
%	H = PROBITNOISE3DPLOT(NOISE, PLOTTYPE, X, Y, MU, VARSIGMA, ...)
%	draws a 3D or contour plot for the probit based classification noise
%	model.
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
%	PROBITNOISEPARAMINIT, NOISE3DPLOT, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



fhandle = str2func(plotType);
  
CZ = cumGaussian((CZ + noise.bias)./sqrt(CZVar));
if nargout > 0
  h = fhandle(CX, CY, CZ, varargin{:});
else
  fhandle(CX, CY, CZ, varargin{:});
end
