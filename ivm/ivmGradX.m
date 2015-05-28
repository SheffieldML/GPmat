function g = ivmGradX(model, x, y);

% IVMGRADX Returns the gradient of the log-likelihood wrt x.
% FORMAT
% DESC returns the gradient of the approximate log likelihood with
% respect to an input location x. This is used for optimising with
% respect to x in the GP-LVM.
% ARG model : the model for which the gradient is being computed.
% ARG x : the input location where the gradient is to be evaluated.
% ARG y : the target location where the gradient is being
% evaluated.
% 
% SEEALSO ivmPosteriorMeanVar, ivmPosteriorGradMeanVar
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005
 
% IVM

[mu, varsigma] = ivmPosteriorMeanVar(model, x);
[dmu, dvs] = ivmPosteriorGradMeanVar(model, x);

g = noiseGradX(model.noise, mu, varsigma, dmu, dvs, y);
