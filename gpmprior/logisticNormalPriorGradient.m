function g = logisticNormalPriorGradient(prior, x)

% LOGISTICNORMALPRIORGRADIENT Gradient wrt x of the logistic-normal prior.

% PRIOR

% Compute gradient of prior
x0 = (x - prior.a) / (prior.b - prior.a);
y = invSigmoid(x0);
g_trans = (x - prior.a) .* (prior.b - x) / (prior.b - prior.a);

g = -(y - prior.mu) / prior.sd^2 ./ (x0 .* (1-x0)) / (prior.b - prior.a) + ...
    -1 ./ g_trans .* (prior.a + prior.b - 2*x) / (prior.b - prior.a);

% invSigmoid(x) = log(x / (1-x))
% invSigmoid'(x) = 1 / (x .* (1-x))
