function l = logisticNormalPriorLogProb(prior, x)

% LOGISTICNORMALPRIORLOGPROB Log probability of logistic-normal prior.

% PRIOR

% Compute log prior

y = invSigmoid((x - prior.a) / (prior.b - prior.a));
%g = sigmoid(y) .* sigmoid(-y) / (prior.b - prior.a);
g = (x - prior.a) .* (prior.b - x) / (prior.b - prior.a)^3;

l = -.5*sum(sum(y.^2/prior.sd^2 + log(2*pi) + ...
		log(prior.sd^2) + 2*log(abs(g))));
