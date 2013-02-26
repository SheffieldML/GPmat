function h = orderedNoisePointPlot(noise, X, y, ...
                              fontName, fontSize, ...
                              markerSize, lineWidth);

% ORDEREDNOISEPOINTPLOT Plot the data-points for the ORDERED noise model.
%
%	Description:
%
%	ORDEREDNOISEPOINTPLOT(NOISE, X, Y, FONTNAME, FONTSIZE, MARKERSIZE,
%	LINEWIDTH) plots the data point locations for the ordered
%	categorical noise structure.
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
%	ORDEREDNOISEPARAMINIT, NOISEPOINTPLOT


%	Copyright (c) 2004, 2005 Neil D. Lawrence



h = [];
symbol = getSymbols(noise.C);
for i = 1:noise.C
  h = [h; plot(X(find(y==i-1), 1), X(find(y==i-1), 2), symbol{i}, 'erasemode', ...
       'xor', 'markerSize', markerSize, 'linewidth', lineWidth)];
  hold on
end

minVals = min([X]);
maxVals = max([X]);

spans = maxVals - minVals;
gaps = spans*.05;

prop = {'xlim', 'ylim'};
for i = 1:2
  set(gca, prop{i}, [minVals(i)-gaps(i) maxVals(i)+gaps(i)]);
end
hold on