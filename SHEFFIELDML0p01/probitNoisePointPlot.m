function h = probitNoisePointPlot(noise, X, y, ...
                              fontName, fontSize, ...
                              markerSize, lineWidth);

% PROBITNOISEPOINTPLOT Plot the data-points for the PROBIT noise model.
%
%	Description:
%
%	PROBITNOISEPOINTPLOT(NOISE, X, Y, FONTNAME, FONTSIZE, MARKERSIZE,
%	LINEWIDTH) plots the data point locations for the probit based
%	classification noise structure.
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
%	PROBITNOISEPARAMINIT, NOISEPOINTPLOT


%	Copyright (c) 2004, 2005 Neil D. Lawrence



pointsNeg = plot(X(find(y(:, 1)==-1), 1), ...
		 X(find(y(:, 1)==-1), 2), ...
		 'gx', 'erasemode', 'xor', ...
		 'markersize', markerSize+2, ...
		 'linewidth', lineWidth);
hold on
pointsPos = plot(X(find(y(:, 1)==1), 1), ...
		 X(find(y(:, 1)==1), 2), 'ro', ...
		 'erasemode', 'xor', ...
		 'markersize', markerSize, ...
		 'linewidth', lineWidth);

h = [pointsNeg; pointsPos];