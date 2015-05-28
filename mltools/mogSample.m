function x = mogSample(model, numSamples)
  
% MOGSAMPLE Sample from a mixture of Gaussians model.
% FORMAT
% DESC samples from a mixture of Gaussians.
% ARG model  : the model that you want to sample from.
% ARG numSamples : the number of samples required.
% RETURN x : the samples from the model.
% 
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : mogCreate

% MLTOOLS
  
p = rand(numSamples, 1);
bins = cumsum(model.prior);
compNo = sum(repmat(p, 1, model.m)<repmat(bins, numSamples, 1), 2);
x = zeros(numSamples, model.d);

for i = 1:model.m
  ind = find(compNo == i);
  x(ind, :) = repmat(model.mean(i, :), length(ind), 1);
  switch model.covtype
   case 'ppca'
    samps = randn(length(ind), model.q);
    samps = samps*model.W{i}';
    samps = samps + randn(length(ind), model.d)*sqrt(model.sigma2(i));
    x(ind, :) = x(ind, :) + samps;
   case 'spherical'
    x(ind, :) = x(ind, :) + randn(length(ind), model.d)* ...
        sqrt(model.sigma(i));
  end
end
