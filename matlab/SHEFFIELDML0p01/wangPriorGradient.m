function g = wangPriorGradient(prior, x)

% WANGPRIORGRADIENT Gradient wrt x of the Wang prior.
%
%	Description:
%	g = wangPriorGradient(prior, x)
%

% Compute gradient of prior
g = -prior.M./x;
