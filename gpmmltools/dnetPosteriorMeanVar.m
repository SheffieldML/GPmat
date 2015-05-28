function [mu, varsigma] = dnetPosteriorMeanVar(model, X);

% DNETPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% FORMAT
% DESC returns the posterior mean and variance for a given set of
% points.
% ARG model : the model for which the posterior will be computed.
% ARG x : the input positions for which the posterior will be
% computed.
% RETURN mu : the mean of the posterior distribution.
% RETURN sigma : the variances of the posterior distributions.
%
% SEEALSO : dnetCreate
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

mu = dnetOut(model, X);

% Return magnification factors instead of variance.
g = modelOutputGradX(model.mapping, X);
varsigma = zeros(size(X, 1), 1);
for n = 1:size(X, 1)
  gTemp = squeeze(g(n, :, :)); 
  varsigma(n) = exp(0.5*logdet(gTemp*gTemp'));
end

%varsigma = repmat(1/model.beta, size(X, 1), 1);
