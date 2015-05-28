function g = pskernelGradient(lntheta, models, prior)

% PSKERNELGRADIENT Gradient on likelihood approximation for point set IVM.

% KERN

% KERN

% PSIVM

if nargin < 3
  prior = 1;
end

g = zeros(size(lntheta));
numTasks = length(models.task);
for taskNo = 1:numTasks
  models.task(taskNo).lntheta = models.lntheta;
  g = g + kernelGradient(lntheta, models.task(taskNo), prior);
end
  
