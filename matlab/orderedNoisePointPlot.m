function orderedNoisePointPlot(noise, X, y, ...
                              fontName, fontSize, ...
                              markerSize, lineWidth);

% ORDEREDNOISEPOINTPLOT Plot the data-points for ordered categorical noise model.

% NOISE

labelColour = {'b', 'g', 'r', 'c', 'm', 'y'};
labelShape = {'o', 'x', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p'};
counter = 0;
while counter < noise.C
  labels{counter+1} = [labelColour{rem(counter, length(labelColour))+1} ...
                     labelShape{rem(counter, length(labelShape))+1}];
  counter = counter +1;
end

numLabels = length(labels);
for i = 1:noise.C
  plot(X(find(y==i-1), 1), X(find(y==i-1), 2), labels{i}, 'erasemode', ...
       'xor', 'markerSize', markerSize, 'linewidth', lineWidth);
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