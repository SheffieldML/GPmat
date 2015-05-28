function varsigma = fgplvmPosteriorVar(model, X);

% FGPLVMPOSTERIORVAR Variances of the posterior at points given by X.
% FORMAT
% DESC returns the posterior mean and variance for a given set of
% points.
% ARG model : the model for which the posterior will be computed.
% ARG x : the input positions for which the posterior will be
% computed.
% RETURN sigma : the variances of the posterior distributions.
%
% SEEALSO : gpPosteriorVar, fgplvmCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009

% FGPLVM

varsigma = gpPosteriorVar(model, X);
