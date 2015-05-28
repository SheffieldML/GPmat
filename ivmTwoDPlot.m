function ivmTwoDPlot(model, iter, X, y)

% IVMTWODPLOT Make a 2-D plot of the IVM.

if nargin < 4
  X = model.X;
  y = model.y;
end
clf
markerSize = 10;
lineWidth = 2;
title(['Iteration ' num2str(iter)])
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

pointUn = plot(X(find(isnan(y(:, 1))), 1), ...
	       X(find(isnan(y(:, 1))), 2), 'm.', ...
	       'erasemode', 'xor', 'markersize', 8);
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 20);

if length(model.I)>0
  % It's probably learnt something.
  xlim = get(gca, 'xlim');
  ylim = get(gca, 'ylim');
  [CX, CY, CZ, CZVar] = ivmMeshVals(model, xlim, ylim, 30);
  noise3dPlot(model.noise, 'ivmContour', CX, CY, CZ, CZVar, lineWidth);   
end
drawnow