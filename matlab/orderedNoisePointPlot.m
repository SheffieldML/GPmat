function orderedNoisePointPlot(noise, X, y, ...
                              fontName, fontSize, ...
                              markerSize, lineWidth);

% ORDEREDNOISEPOINTPLOT Plot the data-points for ordered categorical noise model.

labels{1} = 'bx';
labels{2} = 'ro';
labels{3} = 'g+';
labels{4} = 'ys';
labels{5} = 'mv';
labels{6} = 'c>';
labels{7} = 'w<';

numLabels = length(labels);
for i = 1:noise.C
  plot(X(find(y==i-1), 1), X(find(y==i-1), 2), labels{i}, 'erasemode', ...
       'xor', 'markerSize', markerSize, 'linewidth', lineWidth);
  hold on
end
