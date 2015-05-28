function ncnmRunDataSet3(dataSetName, experimentNo, kernelType, noiseModel, ...
                             selectionCriterion, dVal, seedVal, priorParam);
                             
% NCNMRUNDATASET3 Try the IVM on a data set and save the results.

strPriorParam = num2str(priorParam);
ind = find(strPriorParam == '.');
strPriorParam(ind) = 'p';
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
  [X, y, XTest, yTest] = ncnmLoadData(dataSetName, []);
else
  [X, y, XTest, yTest] = ncnmLoadData(dataSetName, seedVal);
end
labelledIndex = find(~isnan(y));
labelledX = X(labelledIndex, :);
labelledY = y(labelledIndex, :);
capitalName = dataSetName;
capitalName(1) = upper(capitalName(1));
ind = find(capitalName == '.');
capitalName(ind) = 'p';

options = ncnmOptions;
options.extIters = 0;
mu = zeros(size(yTest));
varSigma = zeros(size(yTest));
capName = strtok(capitalName, '_');
load(['../../ivm/matlab/dem' capName '_1_' num2str(seedVal)]);

tic
model = ivmRun(X, y, kernStore, ...
               noiseModel, selectionCriterion, dVal, ...
               options);
  
runTime = toc;
if ~isempty(XTest) & ~isempty(yTest);
  yPred = ivmOut(model, XTest);
  testError = 1-sum(yPred==yTest)/size(yTest, 1);
  fprintf('Labelled data set %s, test error %2.4f\n', dataSetName, testError);
else
  testError = [];
end
[kern, noise, ivmInfo] = ivmDeconstruct(model);
save(['dem' capitalName '_' num2str(experimentNo) ...
      '_' strPriorParam '_' num2str(seedVal)], ...
     'testError', 'ivmInfo', 'kern', 'noise', ...
     'runTime', 'jobId', 'jobName');
