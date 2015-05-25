% ICMTOYREGRESSION Time the point-set IVM and simple sub-sampling.

%/~
importTool('ivm');
%~/

generateMtRegressionData
numTrials = 10;
initTheta = [10 10 10 10];
kernelType = 'rbf';
noiseType = 'gaussian';
selectionCriterion = 'entropy';
display = 1;
innerIters = 50;
outerIters = 5;
prior = 0;

% Run the IVM
dVec = []; %[10 20 30]; %[200 250 300 400 500 600 700 800 900 1000]; 
for dnum = 1:length(dVec)
  d = dVec(dnum);
  fprintf('Size of active set %d\n', d)
  for trials = 1:numTrials
    initTime = cputime;
    models = mtivm(X, y, kernelType, noiseType, 'none', d);
    models.lntheta = log(initTheta);
    models = mtivmOptimise(models, prior, display, innerIters, outerIters);
    paramsIVM{dnum, trials} = kernExtractParam(models.task(1).kern);
    timeIVM(dnum, trials) = cputime - initTime;

    testModels = mtivm(testX, testY, kernelType, noiseType, 'none', d);
    llIVM(dnum, trials) = -mtkernelObjective(paramsIVM{dnum, trials}, models, prior); %testX, ...
%					  testY, 0);
    fprintf('Trial iteration %d complete\n', trials)
    fprintf('Theta %2.4f \n', paramsIVM{dnum, trials})
    fprintf('Likelihood %2.4f Time %2.4f\n', llIVM(dnum, trials), timeIVM(dnum, trials))
  end
end

% Simply sub-sample
iters = 200;
sampsVec = [10 20 30];%[150 200 250 300 350 400 450 500 550 600];
numTasks = length(X);
for sampNum = 1:length(sampsVec); 
  samps = sampsVec(sampNum);
  fprintf('Number of samples %d\n', samps)
  for trials = 1:numTrials
    for task = 1:numTasks
      indices = randperm(size(X{task}, 1));
      Xsamp{task} = X{task}(indices(1:samps), :);
      ysamp{task} = y{task}(indices(1:samps), :);
    end
    initTime = cputime;
    models =  gpPsRun(Xsamp, ysamp, kernelType, noiseType, ...
					initTheta, prior, display, iters);
    paramsSub{sampNum, trials} = kernExtractParam(models.task(1).kern);
    testModels = mtivm(testX, testY, kernelType, noiseType, 'none', []);
    testModels = mtivmInit(testModels);
    llSub(sampNum, trials) = -mtkernelObjective(paramsSub{sampNum, trials}, testModels, prior); %testX, ...

    timeSub(sampNum, trials) = cputime - initTime;
    fprintf('Trial iteration %d complete\n', trials)
    fprintf('Theta %2.4f \n', paramsSub{sampNum, trials})
    fprintf('Likelihood %2.4f Time %2.4f\n', llSub(sampNum, trials), timeSub(sampNum, trials))
    %    displayTasks(Xsamp, ysamp, paramsSub);
  end
end


save('icmlMtRegressionResults.mat', 'timeIVM', 'timeSub', 'llIVM', 'llSub', 'paramsIVM', ...
     'paramsSub', 'dVec', 'sampsVec')


icmlRegressionResults