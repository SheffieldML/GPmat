% DEMUSPSIVM3 Try the ARD IVM on some digits data.

% IVM

dataSetName = 'usps';
experimentNo = 3;

randn('seed', 1e5)
rand('seed', 1e5)

[X, y, XTest, yTest] = mapLoadData(dataSetName);

capitalName = dataSetName;
capitalName(1) = upper(capitalName(1));

options = ivmOptions;
options.kern = {'rbfard', 'lin', 'bias', 'white'};
options.numActive = 500;

mu = zeros(size(yTest));
varSigma = zeros(size(yTest));

a = [ones(4) repmat(2, 4, 4) ...
     repmat(3, 4, 4) repmat(4, 4, 4); ...
     repmat(5, 4, 4) repmat(6, 4, 4) ...
     repmat(7, 4, 4) repmat(8, 4, 4); ...
     repmat(9, 4, 4) repmat(10, 4, 4) ...
     repmat(11, 4, 4) repmat(12, 4, 4); ...
     repmat(13, 4, 4) repmat(14, 4, 4) ...
     repmat(15, 4, 4) repmat(16, 4, 4)];
for i = 1:16
  ardPos = find(a == i)';
  tieIndices{i} = [2+ardPos];
end

tic
% Learn an IVM for each digit
for trainData = 0:9
  index = trainData+1;
  
  % Train the IVM.
  model = ivmRun(X, y(:, index), options);
  
  % Make prediction for this digit.
  [mu(:, index), varSigma(:, index)] = ivmPosteriorMeanVar(model, XTest);
  mu(:, index) = mu(:, index) + model.noise.bias;
  yPred = sign(mu(:, index));
  testError(index) = 1-sum(yPred==yTest(:, index))/size(yTest, 1);
  fprintf('Digit %d, test error %2.4f\n', trainData, testError(index));

  % Deconstruct IVM for saving.
  [kernStore{index}, noiseStore{index}, ...
   ivmInfoStore{index}] = ivmDeconstruct(model);
  save(['dem' capitalName num2str(experimentNo)], 'testError', ...
       'ivmInfoStore', 'kernStore', 'noiseStore')
end
overallTime = toc;

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
save(['dem' capitalName num2str(experimentNo)], 'testError', ...
     'ivmInfoStore', 'kernStore', ...
     'noiseStore', 'overallError', ...
     'confusMat', 'overallTime');
