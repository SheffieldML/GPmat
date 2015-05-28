function [mu, varsigma] = ppcaPosteriorMeanVar(model, X);

% PPCAPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% FORMAT
% DESC returns the posterior mean and variance for a given set of
% points.
% ARG model : the model for which the posterior will be computed.
% ARG x : the input positions for which the posterior will be
% computed.
% RETURN mu : the mean of the posterior distribution.
% RETURN sigma : the variances of the posterior distributions.
%
% SEEALSO : ppcaCreate
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

[mu, phi] = ppcaOut(model, X);
%sqrt(det(phi*model.A*model.A'*phi))
% Include magnification factors here ... need derivatives of outputs ...
% modelOutputGradX ...
varsigma = repmat(1/model.beta, size(X, 1), 1);
