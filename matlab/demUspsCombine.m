% DEMUSPSCOMBINE Select between MLP and RBF and RBF ard network predictors.

importTool('prior');
importTool('kern');
importTool('noise');
importTool('optimi');

[X, y, XTest, yTest] = ivmLoadData('usps');

for j = 1:3
  load(['demUsps' num2str(j)])
  for i = 1:10
    model(j, i) = ivmReconstruct(kernStore{i}, noiseStore{i}, ivmInfoStore{i}, X, y(:, i));
    ll(j, i) = ivmLogLikelihood(model(j, i));
  end
end

[void, index] = min(ll);
mu = zeros(size(XTest, 1), 10);
varSigma = zeros(size(XTest, 1), 10);
for i = 1:10
  [mu(:, i), varSigma(:, i)] = ivmPosteriorMeanVar(model(index(i), i), XTest);
  mu(:, i) = mu(:, i) + model(index(i), i).noise.bias;
end

% Make prediction for all digits.
[void, yPred] = max(mu, [], 2);
[void, yTest] = max(yTest, [], 2);
yPred = yPred - 1;
yTest = yTest - 1;
overallError = 1 - sum(yPred == yTest)/size(yTest, 1);
confusMat = zeros(10);
for i = 1:length(yPred)
  confusMat(yPred(i)+1, yTest(i)+1) = confusMat(yPred(i)+1, yTest(i)+1) + 1;
end
