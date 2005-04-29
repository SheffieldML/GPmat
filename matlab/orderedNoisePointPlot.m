function h = orderedNoisePointPlot(noise, X, y, ...
                              fontName, fontSize, ...
                              markerSize, lineWidth);

% ORDEREDNOISEPOINTPLOT Plot the data-points for ordered categorical noise model.

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