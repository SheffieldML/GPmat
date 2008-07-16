% DEMROBOTWIRELESSNAVIGATE Take some test data for the robot and navigate with it.

% FGPLVM

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


% Initialise X from most likely training point.
X = zeros(size(Y, 1), model.q);
for i = 1:size(Y, 1)
  llVec = fgplvmPointLogLikelihood(model, model.X, repmat(Y(i, :), model.N, 1));
  [void, ind] = max(llVec);
  X(i, :) = model.X(ind, :);
end

% Optimise X
for i =1:size(X, 1)
  Xc(i, :) = fgplvmOptimisePoint(model, X(i, :), Y(i, :), 0, 100);
end

lvmScatterPlot(model, []);
hold on, plot(X(:, 1), X(:, 2), 'bo')
plot(Xc(:, 1), Xc(:, 2), 'gs')

ll(experimentNo) = sum(fgplvmPointLogLikelihood(model, Xc, Y));
fprintf('Test Log-likelihood: %2.4f\n', ll(experimentNo));


