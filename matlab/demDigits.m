% DEMDIGITS Try the IVM on some digits data.

% IVM

load ../data/usps_train
X = ALL_DATA;
yTrain = ALL_T;
dVal = 500;
prior = 0;
display = 0;
kernelType = {'rbfard', 'lin', 'bias', 'white'}
%kernelType = {'rbf', 'bias', 'white'}


noiseModel = 'heaviside';
selectionCriterion = 'entropy';
load ../data/usps_test
xTest = ALL_DATA;
yTest = ALL_T;
u = zeros(size(ALL_T, 1), 10);
varSigma = zeros(size(ALL_T, 1), 10);

a = [ones(4) repmat(2, 4, 4) repmat(3, 4, 4) repmat(4, 4, 4)
     repmat(5, 4, 4) repmat(6, 4, 4) repmat(7, 4, 4) repmat(8, 4, 4)
     repmat(9, 4, 4) repmat(10, 4, 4) ...
     repmat(11, 4, 4) repmat(12, 4, 4)
     repmat(13, 4, 4) repmat(14, 4, 4) ...
     repmat(15, 4, 4) repmat(16, 4, 4)];
for i = 1:16
  ardPos = find(a == i)';
  tieIndices{i} = [2+ardPos];
end

tic
for trainData = 0:9
  y = (yTrain == trainData)*2 - 1;
  model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);
  model.kern = cmpndTieParameters(model.kern, tieIndices);
  model = ivmOptimise(model, prior, display, 100, 4);
  model = ivmOptimiseIVM(model, display);

  
  y = (yTest == trainData)*2 - 1;
  [mu, varSigma(:, trainData+1)] = ivmPosteriorMeanVar(xTest, model);
  u(:, trainData+1) = mu + model.noise.bias;
  yPred = sign(u(:, trainData+1));

  testError(trainData+1) = 1-sum(yPred==y)/size(yTest, 1);
  fprintf('Digit %d, test error %2.4f\n', trainData, testError(trainData+1));
  Istore{trainData+1} = model.I;
  kernStore{trainData+1} = model.kern;
  kernStore{trainData+1} = rmfield(kernStore{trainData+1}, 'Kstore');
  kernStore{trainData+1} = rmfield(kernStore{trainData+1}, 'diagK');
  noiseStore{trainData+1} = model.noise;
  save digitsResults testError Istore kernStore noiseStore
end
overallTime = toc;
[void, yPred] = max(u, [], 2);
yPred = yPred -1;
overallError = 1 - sum(yPred == yTest)/size(yTest, 1);

confusMat = zeros(10);
for i = 1:length(yPred)
  confusMat(yPred(i)+1, ALL_T(i)+1) = confusMat(yPred(i)+1, ALL_T(i)+1) + 1;
end
save digitsResults testError Istore kernStore noiseStore overallError confusMat overallTime
