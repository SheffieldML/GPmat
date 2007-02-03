% DEMCMU35GPLVMRECONSTRUCT2 Reconstruct right leg of CMU 35.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'cmu35gplvm';
experimentNo = 1;

% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
origBias = mean(Y);
origScale = 1./sqrt(var(Y));
%scale = ones(size(scale));
Y = Y - repmat(origBias, size(Y, 1), 1);
Ytest = Ytest - repmat(origBias, size(Ytest, 1), 1);
Y = Y.*repmat(origScale, size(Y, 1), 1);
Ytest = Ytest.*repmat(origScale, size(Ytest, 1), 1);
origYtest = Ytest;
startInd = 63;
Ytest = Ytest(startInd:end, :);
% Indices associated with right leg.
legInd = [8:14];
leglessInd = [1:7 15:size(Ytest, 2)];
YtrainLeg = Y(:, leglessInd);
YtestLeg = Ytest(:, leglessInd);

% Load saved model.
capName = dataSetName;;
capName(1) = upper(capName(1));
load(['dem' capName num2str(experimentNo) '.mat']);

YtrueTest = Ytest;
Ytest(:, legInd) = NaN;
model = gpComputeAlpha(model);
model.dynamics = gpComputeAlpha(model.dynamics);
[Xinit, llo] = fgplvmViterbiSequence(model, model.X_u, Ytest);

Xtest = fgplvmOptimiseSequence(model, Xinit, Ytest, 1, 1000);
Ypred = gpOut(model, Xtest);
