function [mu, varsigma] = fgplvmPosteriorMeanVar(model, X);

% FGPLVMPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.

% FGPLVM

[mu, varsigma] = gpPosteriorMeanVar(model, X);