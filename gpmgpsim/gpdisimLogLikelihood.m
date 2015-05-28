function ll = gpdisimLogLikelihood(model)

% GPDISIMLOGLIKELIHOOD Compute the log likelihood of a GPDISIM model.
% FORMAT
% DESC computes the log likelihood of the given Gaussian process
% for use in a single input motif protein network.
% ARG model : the model for which the log likelihood is computed.
% RETURN ll : the log likelihood of the data set.
% 
% SEEALSO : gpsimCreate, gpsimLogLikeGradient, gpsimObjective
%
% COPYRIGHT : Neil D. Lawrence, 2006

% SHEFFIELDML


dim = size(model.y, 1);
ll = -dim*log(2*pi) - model.logDetK - model.m'*model.invK*model.m;
ll = ll*0.5;

% In case we need priors in.
if isfield(model, 'bprior'),
  ll = ll + kernPriorLogProb(model.kern);
  ll = ll + priorLogProb(model.bprior, model.B);
end

