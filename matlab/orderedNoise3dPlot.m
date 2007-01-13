function h = orderedNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)


% ORDEREDNOISE3DPLOT Draws a 3D or contour plot for the ORDERED noise model.
% FORMAT
% DESC draws a 3D or contour plot for the ordered categorical noise model.
% ARG noise : the noise structure for which the plot is required.
% ARG plotType : string containing the name of the plotting function (for example mesh, contour).
% ARG X : the input X data in the form of a 'mesh' matrix.
% ARG Y : the input Y data in the form of a 'mesh' matrix.
% ARG mu : the input mean in the form of a 'mesh' matrix.
% ARG varSigma : the input variance in the form of a 'mesh' matrix. 
% ARG P1, P2, P3 ... : optional additional arguments for the given plot type.
% RETURN h : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : orderedNoiseParamInit, noise3dPlot, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


CZ = (CZ+noise.bias)./sqrt(CZVar);
fhandle = str2func(plotType);
fhandle(CX, CY, CZ, varargin{:});
hold on
h = [];
for i = 2:noise.C-1
  CZ = CZ - noise.widths(i-1)./sqrt(CZVar);
  h = [h; fhandle(CX, CY, CZ, varargin{:})];
end
