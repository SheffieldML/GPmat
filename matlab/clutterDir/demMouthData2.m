% DEMMOUTHDATA2 Try on Ismael's vowels number 2.

lipsScriptLoadBlabla;

Xtemp = acVec;
scales = sqrt(var(mouthPoints));
for i=1:size(mouthPoints, 2)
  mouthPoints(:, i) = mouthPoints(:, i)/scales(i);
end
ytemp = mouthPoints;

X = Xtemp;
y = ytemp;
%XTest = Xtemp(2:2:end, :);
%yTest = ytemp(2:2:end, :);
dataSetName = 'mouthData';
experimentNo = 2;

randn('seed', 1e5)
rand('seed', 1e5)

capitalName = dataSetName;
capitalName(1) = upper(capitalName(1));

options = ivmOptions;
options.display=1;
options.extIters = 1;
options.noiseIters = 0;
options.kernIters=1000;
kernelType = {'rbfard', 'white'};
noiseType = 'gaussian';
selectionCriterion = 'entropy';


model = ivm(X, y, kernelType, noiseType, selectionCriterion, 500);
model.noise.bias=zeros(size(model.noise.bias));
model = ivmOptimise(model, options);
% Select data-points in an IVM with those kernel parameters.
model = ivmOptimiseIVM(model, options.display);
% Make prediction for the test data.
%[mu, varSigma] = ivmPosteriorMeanVar(model, XTest);
%yPred = mu + repmat(model.noise.bias, [size(mu, 1) 1]);
%testError= sum((yPred-yTest).^2);
%fprintf('Test error %2.4f, seed %2.4f, d %d\n', testError)
      
% Deconstruct IVM for saving.
[kern, noise, ivmInfo] = ivmDeconstruct(model);
save(['dem' capitalName num2str(experimentNo)], ...
         'ivmInfo', 'kern', 'noise')
 