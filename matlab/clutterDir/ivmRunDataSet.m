function ivmRunDataSet(dataSetName, experimentNo, kernelType, noiseModel, ...
                             selectionCriterion, dVal, seedVal);
                             
% IVMRUNDATASET Try the IVM on a data set and save the results.

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

tic
options.kern = kernelType
model = ivmRun(X, y, kernelType, ...
               noiseModel, selectionCriterion, dVal, ...
               options);
  
runTime = toc;
if ~isempty(XTest) & ~isempty(yTest);
  yPred = ivmOut(model, XTest);
  testError = 1-sum(yPred==yTest)/size(yTest, 1);
  fprintf('Data set %s, test error %2.4f\n', dataSetName, testError);
end
[kern, noise, ivmInfo] = ivmDeconstruct(model);
save(['dem' capitalName '_' num2str(experimentNo) '_' num2str(seedVal)], 'testError', ...
     'ivmInfo', 'kern', 'noise', 'runTime', 'jobId', ...
     'jobName')
