function g = wangPriorGradient(prior, x)

% WANGPRIORGRADIENT Gradient wrt x of the Wang prior.

% PRIOR

% Compute gradient of prior
g = -prior.M./x;
