% DEMDIGITS Try the IVM on some digits data.


load ../data/usps_train
X = ALL_DATA;
yTrain = ALL_T;
dVal = 500;
prior = 0;
display = 0;
kernelType = 'rbf';


noiseModel = 'probit';
selectionCriterion = 'entropy';
load ../data/usps_test
xTest = ALL_DATA;
yTest = ALL_T;

u = zeros(size(ALL_T, 1), 10);
varSigma = zeros(size(ALL_T, 1), 10);

tic
for trainData = 0:9
  y = (yTrain == trainData)*2 - 1;
  model = ivmRun(X, y, kernelType, noiseModel, selectionCriterion, dVal, ...
              prior, display, 100, 8);
  
  [mu, varSigma(:, trainData+1)] = ivmPosteriorMeanVar(model, xTest);
  u(:, trainData+1) = mu + model.noise.bias;
  yPred = sign(u(:, trainData+1));
  y = (yTest == trainData)*2 - 1;
  testError(trainData+1) = 1-sum(yPred==y)/size(y, 1);
  fprintf('Digit %d, test error %2.4f\n', trainData, testError(trainData+1));
  [kernStore{trainData+1}, noiseStore{trainData+1}, ivmInfoStore{trainData+1}] = ivmDeconstruct(model);
  save demDigits testError ivmInfoStore kernStore noiseStore
end
overallTime = toc;
[void, yPred] = max(u, [], 2);
yPred = yPred -1;
overallError = 1 - sum(yPred == yTest)/size(yTest, 1);

confusMat = zeros(10);
for i = 1:length(yPred)
  confusMat(yPred(i)+1, ALL_T(i)+1) = confusMat(yPred(i)+1, ALL_T(i)+1) + 1;
end
save demDigits testError ivmInfoStore kernStore noiseStore overallError confusMat overallTime
