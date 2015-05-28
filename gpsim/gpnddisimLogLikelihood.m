function ll = gpnddisimLogLikelihood(model, noprior)

% GPNDDISIMLOGLIKELIHOOD Compute the log likelihood of a GPNDDISIM model.
% FORMAT
% DESC computes the log likelihood of the given Gaussian process
% for use in a single input motif protein network.
% ARG model : the model for which the log likelihood is computed.
% ARG noprior : if true, ignore prior (default: false)
% RETURN ll : the log likelihood of the data set.
% 
% SEEALSO : gpsimCreate, gpsimLogLikeGradient, gpsimObjective
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Jaakko Peltonen, 2011
%
% COPYRIGHT : Antti Honkela, 2012

% GPSIM

if nargin < 2,
  noprior = 0;
end

dim = size(model.y, 1);
ll = -dim*log(2*pi) - model.logDetK - model.m'*model.invK*model.m;
ll = ll*0.5;

if noprior,
  return;
end

% In case we need priors in.
if isfield(model, 'bprior'),
  ll = ll + kernPriorLogProb(model.kern);
  if model.numGenes>0,
    ll = ll + priorLogProb(model.bprior, model.B);
    ll = ll + priorLogProb(model.simMeanPrior, model.simMean);
  end;
end

if isfield(model, 'disimStartMeanPrior'),
  ll = ll + priorLogProb(model.disimStartMeanPrior, ...
			 model.disimStartMean);
end
