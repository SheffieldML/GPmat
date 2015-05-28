function h = orderedNoisePointPlot(noise, X, y, ...
                              fontName, fontSize, ...
                              markerSize, lineWidth);

% ORDEREDNOISEPOINTPLOT Plot the data-points for the ORDERED noise model.
% FORMAT
% DESC plots the data point locations for the ordered categorical
%  noise structure.
% ARG noise : the noise structure which is to be plotted.
% ARG X : input locations to be plotted.
% ARG y : target locations to be plotted.
% ARG fontName : name of any fonts to be used.
% ARG fontSize : size of any fonts to be used.
% ARG markerSize : size of the makers to be used.
% ARG lineWidth : size of any lines to be plotted.
%
% SEEALSO : orderedNoiseParamInit, noisePointPlot
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


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