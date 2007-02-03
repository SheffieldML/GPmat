function ll = gpPointLogLikelihood(model, x, y)

% GPPOINTLOGLIKELIHOOD Log-likelihood of a test point for a GP.
% FORMAT
% DESC returns the log likelihood of a latent point and an observed
% data point for the posterior prediction of a GP model.
% ARG model : the model for which the point prediction will be
% made.
% ARG x : the input point for which the posterior distribution
% will be evaluated.
% ARG y : the target point for which the posterior
% distribution will be evaluated.
%
% SEEALSO : gpCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GP

logTwoPi = log(2*pi);
[mu, varSigma] = gpPosteriorMeanVar(model, x);
ll = zeros(size(x, 1), 1);
ydiff = y-mu;
ll = log(varSigma) + (ydiff.*ydiff)./varSigma +logTwoPi;
ll(find(isnan(ll)))=0;
ll = -0.5*sum(ll, 2);
