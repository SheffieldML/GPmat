% ICMTOYREGRESSION Time the point-set IVM and simple sub-sampling.

generatePsRegressionData
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
    models = psivm(X, y, kernelType, noiseType, 'none', d);
    models.lntheta = log(initTheta);
    models = psivmOptimise(models, prior, display, innerIters, outerIters);
    thetaIVM{dnum, trials} = exp(models.lntheta);
    timeIVM(dnum, trials) = cputime - initTime;

    testModels = psivm(testX, testY, kernelType, noiseType, 'none', d);
    llIVM(dnum, trials) = -pskernelObjective(log(thetaIVM{dnum, trials}), models, prior); %testX, ...
%					  testY, 0);
    fprintf('Trial iteration %d complete\n', trials)
    fprintf('Theta %2.4f \n', thetaIVM{dnum, trials})
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
    thetaSub{sampNum, trials} = exp(models.lntheta);
    testModels = psivm(testX, testY, kernelType, noiseType, 'none', []);
    testModels.lntheta = models.lntheta;
    testModels = psivmInit(testModels);
    llSub(sampNum, trials) = -pskernelObjective(log(thetaSub{sampNum, trials}), testModels, prior); %testX, ...

    timeSub(sampNum, trials) = cputime - initTime;
    llSub(sampNum, trials) = -pskernelObjective(log(thetaSub{sampNum, ...
		    trials}), testModels);
    fprintf('Trial iteration %d complete\n', trials)
    fprintf('Theta %2.4f \n', thetaSub{sampNum, trials})
    fprintf('Likelihood %2.4f Time %2.4f\n', llSub(sampNum, trials), timeSub(sampNum, trials))
    %    displayTasks(Xsamp, ysamp, thetaSub);
  end
end


save('icmlToyRegressionResults.mat', 'timeIVM', 'timeSub', 'llIVM', 'llSub', 'thetaIVM', ...
     'thetaSub', 'dVec', 'sampsVec')


icmlRegressionResults