function g = invgammaPriorGradient(prior, x)

% INVGAMMAPRIORGRADIENT Gradient wrt x of the log Gaussian prior.

% IVM

% Compute gradient of prior
D = length(x);
g =-(prior.a+1)./x + prior.b./(x.*x);
