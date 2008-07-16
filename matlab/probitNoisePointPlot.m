function h = probitNoisePointPlot(noise, X, y, ...
                              fontName, fontSize, ...
                              markerSize, lineWidth);

% PROBITNOISEPOINTPLOT Plot the data-points for the PROBIT noise model.
% FORMAT
% DESC plots the data point locations for the probit based classification
%  noise structure.
% ARG noise : the noise structure which is to be plotted.
% ARG X : input locations to be plotted.
% ARG y : target locations to be plotted.
% ARG fontName : name of any fonts to be used.
% ARG fontSize : size of any fonts to be used.
% ARG markerSize : size of the makers to be used.
% ARG lineWidth : size of any lines to be plotted.
%
% SEEALSO : probitNoiseParamInit, noisePointPlot
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


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