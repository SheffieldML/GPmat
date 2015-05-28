function f = pskernelObjective(lntheta, models, prior)

% PSKERNELOBJECTIVE Likelihood approximation for point set IVM.

% KERN

% KERN

% PSIVM

if nargin < 3
  prior = 1;
end
f = 0;
numTasks = length(models.task);
for taskNo = 1:numTasks
  models.task(taskNo).lntheta = models.lntheta;
  f = f + kernelObjective(lntheta, models.task(taskNo), prior);
end
