function [mu, varsigma] = fgplvmPosteriorMeanVar(model, X);

% FGPLVMPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
%
% [mu, varsigma] = fgplvmPosteriorMeanVar(model, X);
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmPosteriorMeanVar.m version 1.1



[mu, varsigma] = gpPosteriorMeanVar(model, X);