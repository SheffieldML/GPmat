function M = gplvmdatamakeavi1d(X, Y, invK, theta, YLbls, meanData, activeSet, ...
			visualiseFunction, visualiseModify, varargin)

% GPLVMDATAMAKEAVI1D Create a movie of the data along the line.

numData = size(Y, 1);
if nargin < 6
  YLbls = ones(numData, 1);
end

[void, order] = sort(X);
Y = Y(order, :);
figure
set(gcf, 'menubar', 'none')
set(gcf, 'position', [232 572 115 115/varargin{1}(1)*varargin{1}(2)])
imageAxesa = subplot(1, 1, 1);
set(imageAxesa, 'position', [0 0 varargin{1}(1)/115  varargin{1}(1)/115]);
visHandle = feval(visualiseFunction, Y(1, :), varargin{:});
set(visHandle, 'erasemode', 'xor')
colorMap gray
set(imageAxesa, 'visible', 'off')
for i = 1:size(Y, 1)
  feval(visualiseModify, visHandle, ...
	Y(i, :)+meanData, varargin{:});
  M(i) = getframe;
end