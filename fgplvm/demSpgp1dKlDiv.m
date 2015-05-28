% DEMSPGP1DKLDIV Compute the KL divergences between the demos and the truth.

% FPGPLVM

printDiagram = 1;
demSpgp1d1
demSpgp1d2
demSpgp1d3
demSpgp1d4
load demSpgp1d4
[truem, truecovar] = gpPosteriorMeanCovar(model, model.X);
tmodel = model;
for i = 1:3
  load(['demSpgp1d' num2str(i)]);
  
  [m, covar] = gpPosteriorMeanCovar(model, model.X);
  KL1(i) = kldivGaussian(truem, truecovar, m, covar);
  KL2(i) = kldivGaussian(m, covar, truem, truecovar);  
end

