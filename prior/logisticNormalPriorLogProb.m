function l = logisticNormalPriorLogProb(prior, x)

% LOGISTICNORMALPRIORLOGPROB Log probability of logistic-normal prior.

% PRIOR

% Compute log prior

y = invSigmoid((x - prior.a) / (prior.b - prior.a));
%g = sigmoid(y) .* sigmoid(-y) / (prior.b - prior.a);
g = (x - prior.a) .* (prior.b - x) / (prior.b - prior.a);

if any(x <= prior.a | x >= prior.b),
  l = -Inf;
else
  l = -.5*sum(sum((y - prior.mu).^2/prior.sd^2 + log(2*pi) + ...
		  log(prior.sd^2) + 2*log(abs(g))));
end
