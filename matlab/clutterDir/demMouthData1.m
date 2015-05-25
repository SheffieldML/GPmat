% DEMMOUTH1 Try on Ismaels vowel's.


dataSetName = 'mouthData';
experimentNo = 1;

randn('seed', 1e5)
rand('seed', 1e5)

capitalName = dataSetName;
capitalName(1) = upper(capitalName(1));

options = ivmOptions;
options.display=1;
options.extIters = 1;
options.kernIters=500;
kernelType = {'mlpard', 'bias', 'white'};
noiseType = 'gaussian';
selectionCriterion = 'entropy';


[X, y, XTest, yTest] = ivmLoadData(dataSetName);
model = ivm(X, y, kernelType, noiseType, selectionCriterion, 200);
model = ivmOptimise(model, options);
% Select data-points in an IVM with those kernel parameters.
model = ivmOptimiseIVM(model, options.display);
% Make prediction for the test data.
[mu, varSigma] = ivmPosteriorMeanVar(model, XTest);
yPred = mu + repmat(model.noise.bias, [size(mu, 1) 1]);
testError= sum((yPred-yTest).^2);
fprintf('Test error %2.4f, seed %2.4f, d %d\n', testError)
      
% Deconstruct IVM for saving.
[kern, noise, ivmInfo] = ivmDeconstruct(model);
save(['dem' capitalName num2str(experimentNo)], 'testError', ...
         'ivmInfo', 'kern', 'noise')
 