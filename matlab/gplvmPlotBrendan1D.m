% GPLVMPLOTBRENDAN1D Generate a figure containing fantasy brendans for NIPS paper.

% load data
load frey_rawface.mat

Y = double(ff)';

meanData = mean(Y);
Y = Y  - repmat(meanData, size(Y, 1), 1);

load gplvmBrendan1D.mat

[K, invK] = computeKernel(X(activeSet, :), theta);

% Visualise the results
numFrames = 64;
activeX = X(activeSet, :);
activeY = Y(activeSet, :);

numData = size(Y, 1);

XTest = linspace(min(X(:, 1)), max(X(:, 1)), numFrames)';
testY = manifoldOutputs(XTest, activeX, activeY, theta, invK);
testY = testY + repmat(meanData, numFrames, 1);
testY = testY - min(min(testY));
testY = testY/max(max(testY));

% Display the uniformly spaced fantasies of Brendan
for i = 1:numFrames
  subplot(4, 16, i)
  image(reshape(testY(i, :)*64, 20, 28)');
  colormap gray
  axis image
  axis off
end
set(gcf, 'paperunits', 'normalized')
set(gcf, 'position', [20 500 800 200]);
set(gcf, 'paperposition', [0.1 0.4 0.8 0.2])
figure
% Now find the closed true brendan's and plot them
for i = 1:numFrames
  [void, closest] = min(dist2(XTest(i, :), X));
  subplot(4, 16, i)
  image(reshape((Y(closest, :)+meanData), 20, 28)'/4);
  colormap gray
  axis image
  axis off
end
set(gcf, 'paperunits', 'normalized')
set(gcf, 'position', [70 450 800 200]);
set(gcf, 'paperposition', [0.1 0.4 0.8 0.2])