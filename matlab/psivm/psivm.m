function models = psivm(X, y, kernelType, noiseType, selectionCriterion, d)

% PSIVM Initialise an point set ivm model.

% PSIVM

if nargin < 6
  d = [];
end
models.kernelType = kernelType;
models.noiseType = noiseType;
models.selectionCriteria = selectionCriterion;

numTasks = length(X);
for taskNo = 1:numTasks
  models.task(taskNo) = ...
      ivm(X{taskNo}, y{taskNo}, kernelType, ...
	  noiseType, selectionCriterion, d);
end
models.lntheta = models.task(1).lntheta;
models.d = d;
% Remove fields from the sub-structures.
for taskNo = 1:numTasks
  rmfield(models.task, 'lntheta');
  rmfield(models.task, 'd');
  rmfield(models.task, 'kernelType');
end
