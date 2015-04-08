function ivmRunDataSetRegression2(dataSetName, experimentNo, kernelType, noiseModel, ...
                             selectionCriterion, dVal, seedVal);
                             
% IVMRUNDATASETREGRESSION2 Try the IVM on a data set and save the results.

if nargin < 7
  seedVal = [];
  randn('seed', 1e5)
  rand('seed', 1e5)
else
  randn('seed', seedVal);
  rand('seed', seedVal);
end

jobId = getenv('PBS_JOBID');
jobName = getenv('PBS_JOBNAME');
fprintf('Seed %2.0e\n', seedVal);

if isempty(seedVal)
  [X, y, XTest, yTest] = mapLoadData(dataSetName);
else
  [X, y, XTest, yTest] = mapLoadData(dataSetName, seedVal);
end

capitalName = dataSetName;
capitalName(1) = upper(capitalName(1));

options = ivmOptions;

mu = zeros(size(yTest));
varSigma = zeros(size(yTest));

kern = kernCreate(X, kernelType);
kern.comp{1}.inputScales = repmat(0.5, size(kern.comp{1}.inputScales));

tic
model = ivmRun(X, y, kern, ...
               noiseModel, selectionCriterion, dVal, ...
               options);
  
runTime = toc;
if ~isempty(XTest) & ~isempty(yTest);
  [mu, varSigma] = ivmPosteriorMeanVar(model, XTest);
  yPred = mu + model.noise.bias;
  testError = 0.5*mean((yPred-yTest).^2);
  fprintf('Test error %2.4f, seed %2.4f, d %d\n', testError, ...
          seedVal, dVal);
end
[kern, noise, ivmInfo] = ivmDeconstruct(model);
save(['dem' capitalName '_' num2str(experimentNo) '_d' num2str(dVal) '_seed' num2str(seedVal)], 'testError', ...
     'ivmInfo', 'kern', 'noise', 'runTime', 'jobId', ...
     'jobName')
