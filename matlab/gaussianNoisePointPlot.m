function h = gaussianNoisePointPlot(noise, X, y, ...
                              fontName, fontSize, ...
                              markerSize, lineWidth);

% GAUSSIANNOISEPOINTPLOT Plot the data-points for ordered categorical noise model.

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