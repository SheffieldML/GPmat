% DEMCMU49BALANCEARM5 Demonstrate latent force model on CMU data with regression.

% MULTIGP

rand('seed', 1e6);
randn('seed', 1e6);

dataSetName = 'cmu49BalanceArm';
experimentNo = 5;


% load data
[y, void, yTest, void] = lvmLoadData(dataSetName);

scaleVal = sqrt(sum(var(y)));
y = y/scaleVal;
yTest = yTest/scaleVal;

xTrain = y(:, 1:3);
yTrain = y(:, 4:end);

xTest = yTest(:, 1:3);
yTest = yTest(:, 4:end);

% Set up the model
options = gpOptions('ftc');

% Use the full Gaussian process model.
q = size(xTrain, 2);
d = size(yTrain, 2);

iters = 1000;
display = 1;

for i = 1:d
  model{i} = gpCreate(q, 1, xTrain, yTrain(:, i), options);
  % Optimise the model.
  model{i} = gpOptimise(model{i}, display, iters);
end

% Save results
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');
