function gaussianNoisePointPlot(noise, X, y, ...
                              fontName, fontSize, ...
                              markerSize, lineWidth);

% GAUSSIANNOISEPOINTPLOT Plot the data-points for ordered categorical noise model.

plot3(X(:, 1), X(:, 2), y, 'b.', 'erasemode', 'xor',  'markerSize', markerSize, 'linewidth', lineWidth);
hold on