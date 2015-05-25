function gplvmvisualiseData1D(X, Y, invK, theta, YLbls, meanData, activeSet, ...
			visualiseFunction, visualiseModify, varargin)

% GPLVMVISUALISEDATA1D Visualise the ordering of the data across the line.

numData = size(Y, 1);
if nargin < 6
  YLbls = ones(numData, 1);
end

[void, order] = sort(X);
Y = Y(order, :);
figure
imageAxesa = subplot(1, 1, 1);
visHandle = feval(visualiseFunction, Y(1, :), varargin{:});
set(visHandle, 'erasemode', 'xor')
colorMap gray
set(imageAxesa, 'visible', 'off')
for i = 1:size(Y, 1)
    feval(visualiseModify, visHandle, ...
	  Y(i, :)+meanData, varargin{:});
    pause(0.1)
end