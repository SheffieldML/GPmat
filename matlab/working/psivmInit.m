function models = psivmInit(models)

% PSIVMINIT Initialise point-set IVM.

% PSIVM

for i = 1:length(models.task)
  models.task(i) = ivmInit(models.task(i));
end
