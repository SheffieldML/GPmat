function [mu, varsigma] = fgplvmPosteriorMeanVar(model, X);

% FGPLVMPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% FORMAT
% DESC returns the posterior mean and variance for a given set of
% points.
% ARG model : the model for which the posterior will be computed.
% ARG x : the input positions for which the posterior will be
% computed.
% RETURN mu : the mean of the posterior distribution.
% RETURN sigma : the variances of the posterior distributions.
%
% SEEALSO : gpPosteriorMeanVar, fgplvmCreate
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% FGPLVM

[mu, varsigma] = gpPosteriorMeanVar(model, X);
