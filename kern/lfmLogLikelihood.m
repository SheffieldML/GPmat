function ll = lfmLogLikelihood(model)

% LFMLOGLIKELIHOOD Compute the log likelihood of a LFM model.
% FORMAT
% DESC computes the log likelihood of the given Gaussian process
% for use in a latent force model.
% ARG model : the model for which the log likelihood is computed.
% RETURN ll : the log likelihood of the data set.
% 
% SEEALSO : lfmCreate, lfmLogLikeGradient, lfmObjective
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

dim = size(model.y, 1);
ll = -dim*log(2*pi) - model.logDetK - model.m'*model.invK*model.m;
ll = ll*0.5;
%/~ In case we need priors in.
%ll = ll + kernPriorLogProb(model.kern);
%ll = ll + priorLogProb(model.bprior, model.B);
%~/
