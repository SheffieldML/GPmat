function h = gaussianNoisePointPlot(noise, X, y, fontName, fontSize, ...
                                    markerSize, lineWidth)

% GAUSSIANNOISEPOINTPLOT Plot the data-points for the GAUSSIAN noise model.
%
%	Description:
%
%	GAUSSIANNOISEPOINTPLOT(NOISE, X, Y, FONTNAME, FONTSIZE, MARKERSIZE,
%	LINEWIDTH) plots the data point locations for the Gaussian noise
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
%	GAUSSIANNOISEPARAMINIT, NOISEPOINTPLOT


%	Copyright (c) 2004, 2005 Neil D. Lawrence



h = plot3(X(:, 1), X(:, 2), y, 'r.', 'erasemode', 'xor',  'markerSize', markerSize, 'linewidth', lineWidth);

minVals = min([X y]);
maxVals = max([X y]);

spans = maxVals - minVals;
gaps = spans*.05;

prop = {'xlim', 'ylim', 'zlim'};
for i = 1:3
  set(gca, prop{i}, [minVals(i)-gaps(i) maxVals(i)+gaps(i)]);
end
hold on