function vivmRunDataSet(dataSetName, invariance, experimentNo, ...
                             selectionCriterion, dVal, seedVal);
                             
% VIVMRUNDATASET Try the virtual IVM on a data set and save the results.

% IVM

if nargin < 6
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
options.extIters = 0; % don't optimise parameters
mu = zeros(size(yTest));
varSigma = zeros(size(yTest));

load(['dem' capitalName '_' ...
      num2str(experimentNo) '_' ...
      num2str(seedVal) '.mat'])

% Create new training set with virtual SVs.
origIndex = ivmInfoStore.I;
[X, y] = ivmVirtual(X(ivmInfoStore.I, :), y(ivmInfoStore.I, ...
                                                  :), invariance);
tic
model = ivmRun(X, y, kernStore, ...
               noiseStore, selectionCriterion, dVal, ...
               options);
  
runTime = toc;
if ~isempty(XTest) & ~isempty(yTest);
  yPred = ivmOut(model, XTest);
  testError = 1-sum(yPred==yTest)/size(yTest, 1);
  fprintf('Data set %s, test error %2.4f\n', dataSetName, testError);
end
[kern, noise, ivmInfo] = ivmDeconstruct(model);
save(['dem' capitalName '_' invariance '_' num2str(experimentNo) '_d' num2str(dVal) '_seed' num2str(seedVal)], 'testError', ...
     'ivmInfo', 'kern', 'noise', 'runTime', 'jobId', ...
     'jobName', 'origIndex')

