function [varsigma] = ppcaPosteriorVar(model, X);

% PPCAPOSTERIORVAR Mean and variances of the posterior at points given by X.
% FORMAT
% DESC returns the posterior mean and variance for a given set of
% points.
% ARG model : the model for which the posterior will be computed.
% ARG x : the input positions for which the posterior will be
% computed.
% RETURN sigma : the variances of the posterior distributions.
%
% SEEALSO : ppcaCreate, ppcaPosteriorMeanVar
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

varsigma = repmat(1/model.beta, size(X, 1), 1);
