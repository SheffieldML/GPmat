function h = noisePointPlot(noise, X, y,  ...
                        fontName, fontSize, ...
                        markerSize, lineWidth)

% NOISEPOINTPLOT Plot the data-points for the given noise model.
%
%	Description:
%
%	NOISEPOINTPLOT(NOISE, X, Y, FONTNAME, FONTSIZE, MARKERSIZE,
%	LINEWIDTH) plots the data point locations for the given noise
%	structure.
%	 Arguments:
%	  NOISE - the noise structure which is to be plotted.
%	  X - input locations to be plotted.
%	  Y - target locations to be plotted.
%	  FONTNAME - name of any fonts to be used.
%	  FONTSIZE - size of any fonts to be used.
%	  MARKERSIZE - size of the makers to be used.
%	  LINEWIDTH - size of any lines to be plotted.
%	
%
%	See also
%	NOISEPARAMINIT, NOISECREATE


%	Copyright (c) 2004, 2005 Neil D. Lawrence


fhandle = str2func([noise.type 'NoisePointPlot']);
h = fhandle(noise, X, y, ...
      fontName, fontSize, ...
      markerSize, lineWidth);
if ~isoctave
  set(gca, 'fontname', fontName, 'fontsize', fontSize);
end