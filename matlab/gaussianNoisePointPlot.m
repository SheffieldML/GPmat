function h = gaussianNoisePointPlot(noise, X, y, fontName, fontSize, ...
                                    markerSize, lineWidth)

% GAUSSIANNOISEPOINTPLOT Plot the data-points for the GAUSSIAN noise model.
% FORMAT
% DESC plots the data point locations for the Gaussian
%  noise structure.
% ARG noise : the noise structure which is to be plotted.
% ARG X : input locations to be plotted.
% ARG y : target locations to be plotted.
% ARG fontName : name of any fonts to be used.
% ARG fontSize : size of any fonts to be used.
% ARG markerSize : size of the makers to be used.
% ARG lineWidth : size of any lines to be plotted.
%
% SEEALSO : gaussianNoiseParamInit, noisePointPlot
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


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