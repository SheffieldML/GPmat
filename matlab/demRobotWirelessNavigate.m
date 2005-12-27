% DEMROBOTWIRELESSNAVIGATE Take some test data for the robot and navigate with it.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'robotWireless';
experimentNo = 2;

% Load results of learning.
capName = dataSetName;;
capName(1) = upper(capName(1));
load(['dem' capName num2str(experimentNo) '.mat']);


% load data
[Y, lbls] = lvmLoadData([dataSetName 'Test']);


X = mlpOut(model.back, Y);

for i =1:size(X, 1)
  Xc(i, :) = fgplvmOptimisePoint(model, X(i, :), Y(i, :), 1, 100);
end
  
% Plot for Dieter. Note that we are not using Xc.
lvmScatterPlot(model);
hold on, plot(X(:, 1), X(:, 2), 'bo')