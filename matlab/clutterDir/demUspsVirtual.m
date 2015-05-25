% DEMUSPSVIRTUAL Use virtual informative vectors to train USPS data.

dataSetName = 'usps';
experimentNo = 1;

randn('seed', 1e5)
rand('seed', 1e5)

[origX, origy, XTest, yTest] = ivmLoadData(dataSetName);
load demUsps1.mat


capitalName = dataSetName;
capitalName(1) = upper(capitalName(1));
dVal = 1000;

selectionCriterion = 'entropy';

mu = zeros(size(yTest));
varSigma = zeros(size(yTest));

clear X
baseVal = min(min(origX));
for trainData = 0:9

  % Create new training set with virtual SVs.
  index = trainData + 1;
  numActive = length(ivmInfoStore{index}.I);
  X{index} = zeros(numActive*5, size(origX, 2));
  y{index} = zeros(numActive*5, size(origy, 2));
  for j = 1:length(ivmInfoStore{index}.I)
    xPoint =  origX(ivmInfoStore{index}.I(j), :);
    X{index}(5*j-4, :) = xPoint;
    Xmat = reshape(xPoint, 16, 16);
    XmatU = [Xmat(2:end, :); repmat(baseVal, 1, 16)];
    X{index}(5*j-3, :) = XmatU(:)';
    XmatD = [repmat(baseVal, 1, 16); Xmat(1:end-1, :)];
    X{index}(5*j-2, :) = XmatD(:)';
    XmatL = [Xmat(:, 2:end) repmat(baseVal, 16, 1)];
    X{index}(5*j-1, :) = XmatL(:)';
    XmatR = [repmat(baseVal, 16, 1) Xmat(:, 1:end-1)];
    X{index}(5*j, :) = XmatR(:)';
    for i = 0:4
      y{index}(5*j-i, :) = origy(ivmInfoStore{index}.I(j), :);
    end
  end
  
  % Train the IVM.
  % Select data-points in an IVM with those kernel parameters.
  model = ivm(X{index}, y{index}(:, index), {'rbf', 'lin', 'bias', 'white'}, noiseStore{index}.type, selectionCriterion, dVal);
  model.kern = kernStore{index};
  model.noise = noiseStore{index};
  ivmDisplay(model);
  model = ivmOptimiseIVM(model, 0);
  
  % Make prediction for this digit.
  [mu(:, index), varSigma(:, index)] = ivmPosteriorMeanVar(model, XTest);
  mu(:, index) = mu(:, index) + model.noise.bias;
  yPred = sign(mu(:, index));
  testError(index) = 1-sum(yPred==yTest(:, index))/size(yTest, 1);
  fprintf('Digit %d, test error %2.4f\n', trainData, testError(index));

  % Deconstruct IVM for saving.
  [kernStore{index}, noiseStore{index}, ...
   ivmInfoStore{index}] = ivmDeconstruct(model);
  save(['dem' capitalName 'Virtual' num2str(experimentNo)], 'testError', ...
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
save(['dem' capitalName 'Virtual' num2str(experimentNo)], 'testError', ...
     'ivmInfoStore', 'kernStore', ...
     'noiseStore', 'overallError', ...
     'confusMat', 'overallTime');
