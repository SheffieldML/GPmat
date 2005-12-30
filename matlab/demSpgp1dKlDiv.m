% DEMSPGP1DKLDIV Compute the KL divergences between the demos and the truth.

% FPGPLVM

load demSpgp1d4
[truem, truecovar] = gpPosteriorMeanCovar(model, model.X);
for i = 1:3
  load(['demSpgp1d' num2str(i)]);
  [m, covar] = gpPosteriorMeanCovar(model, model.X);
  KL1(i) = kldivGaussian(truem, truecovar, m, covar);
  KL2(i) = kldivGaussian(m, covar, truem, truecovar);
  
end

