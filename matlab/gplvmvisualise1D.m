function gplvmvisualise1D(X, Y, invK, theta, YLbls, meanData, activeSet, ...
			visualiseFunction, visualiseModify, varargin)

% GPLVMVISUALISE1D Visualise the fantasies along a line (as a movie).

numFrames = 10000;
activeX = X(activeSet, :);
activeY = Y(activeSet, :);

numData = size(Y, 1);
if nargin < 6
  YLbls = ones(numData, 1);
end

XTest = linspace(min(X(:, 1)), max(X(:, 1)), numFrames)';
[testY, testYVar] = manifoldOutputs(XTest, activeX, activeY, theta, invK);
vars = sort(testYVar);
pickPoint = ceil(0.75*numFrames);
maxVar = vars(pickPoint);
figure(1)
imageAxesa = subplot(1, 1, 1);
visHandle = feval(visualiseFunction, Y(1, :), varargin{:});
set(visHandle, 'erasemode', 'xor')
colorMap gray
set(imageAxesa, 'visible', 'off')
for i = 1:size(testY, 1)
  if testYVar(i) < maxVar
    feval(visualiseModify, visHandle, ...
	  testY(i, :) + meanData, varargin{:});
    pause(0.005)
  end
end