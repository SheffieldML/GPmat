function h = noisePointPlot(noise, X, y,  ...
                        fontName, fontSize, ...
                        markerSize, lineWidth)

% NOISEPOINTPLOT 

% NOISE

fhandle = str2func([noise.type 'NoisePointPlot']);
h = fhandle(noise, X, y, ...
      fontName, fontSize, ...
      markerSize, lineWidth);

set(gca, 'fontname', fontName, 'fontsize', fontSize);