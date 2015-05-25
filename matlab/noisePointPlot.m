function h = noisePointPlot(noise, X, y,  ...
                        fontName, fontSize, ...
                        markerSize, lineWidth)

% NOISEPOINTPLOT Plot the data-points for the given noise model.
% FORMAT
% DESC plots the data point locations for the given
%  noise structure.
% ARG noise : the noise structure which is to be plotted.
% ARG X : input locations to be plotted.
% ARG y : target locations to be plotted.
% ARG fontName : name of any fonts to be used.
% ARG fontSize : size of any fonts to be used.
% ARG markerSize : size of the makers to be used.
% ARG lineWidth : size of any lines to be plotted.
%
% SEEALSO : noiseParamInit, noiseCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE

fhandle = str2func([noise.type 'NoisePointPlot']);
h = fhandle(noise, X, y, ...
      fontName, fontSize, ...
      markerSize, lineWidth);
if ~isoctave
  set(gca, 'fontname', fontName, 'fontsize', fontSize);
end
