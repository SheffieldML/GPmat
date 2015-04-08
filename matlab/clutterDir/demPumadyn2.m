% DEMPUMADYN2 Try the IVM on the pumadyn robot arm data.


dataSetName = 'pumadyn';
experimentNo = 2;

randn('seed', 1e5)
rand('seed', 1e5)

capitalName = dataSetName;
capitalName(1) = upper(capitalName(1));

options = ivmOptions;
options.extIters = 4;
kernelType = {'rbfard', 'bias', 'white'};
noiseModel = 'gaussian';
options.kern = kernelType;
options.noise = noiseModel;


dVals = [50 100 200 500 1000];
seedVals = [8:10]*1e5;
testError = zeros(length(seedVals), length(dVals));
overallTime = testError;
for i = 1:length(seedVals);
  [X, y, XTest, yTest] = mapLoadData(dataSetName, seedVals(i));
  Xscale{i} = sqrt(var(X));
  yscale{i} = sqrt(var(y));
  for j = 1:length(Xscale{i});
    X(:, j) = X(:, j)/Xscale{i}(j);
    XTest(:, j) = XTest(:, j)/Xscale{i}(j);
  end
  y = y(:, 1)/yscale{i};
  
  for j = 1:length(dVals);

    tic;
    % Train the IVM.
    options.numActive = dVals(j);
    model = ivmRun(X, y, options);
    
    overallTime(i, j) = toc;
    % Make prediction for the test data.
    [mu, varSigma] = ivmPosteriorMeanVar(model, XTest);
    yPred = mu + model.noise.bias;
    testError(i, j) = sum((yPred*yscale{i}-yTest).^2);
    fprintf('Test error %2.4f, seed %2.4f, d %d\n', testError(i, j), ...
            seedVals(i), dVals(j));
      
    % Deconstruct IVM for saving.
    [kernStore{i, j}, noiseStore{i, j}, ivmInfoStore{i, j}] = ivmDeconstruct(model);
    save(['dem' capitalName num2str(experimentNo)], 'testError', ...
         'ivmInfoStore', 'kernStore', 'noiseStore')
  end
end