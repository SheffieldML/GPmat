function [gParam, gX_u, gX] = gpsimCandidateLogLikeGradients(model, X, M, X_u)

% GPSIMCANDIDATELOGLIKEGRADIENTS Compute the gradients for the parameters of candidate genes.
% FORMAT
% DESC computes the gradients of the Gaussian process posterior 
% likelihood with respect to the parameters of the model.
% ARG : model : the model structure for which gradients are computed.
% RETURN gParam : the gradient of the log likelihood with respect to
% the model parameters.
%
% SEEALSO : gpsimAddCandidate, gpsimCandidateLogLikelihood, modelLogLikeGradients, fgplvmLogLikeGradients
%
% COPYRIGHT : Neil D. Lawrence, 2007

% SHEFFIELDML

if nargin < 4
  if nargin < 3
    M = model.candidate.m;
  end
  if nargin < 2
    t = model.candidate.t;
  end
end

gX_u = [];
gX = [];

%g_scaleBias = gpScaleBiasGradient(model);
if isfield(model, 'meanFunction') & ~isempty(model.meanFunction)
  g_meanFunc = gpMeanFunctionGradient(model);
else
  g_meanFunc = [];
end

[gK_uf, gK_star] = gpsimCandidateCovGrads(model, M);
  
%%% Compute Gradients of Kernel Parameters %%%
g_param = kernGradient(model.candidate.kern, model.candidate.t, model.t, gK_uf);
  
% deal with diagonal term's affect on kernel parameters.
g_param = g_param + kernGradient(model.candidate.kern, model.t, gK_star);

% Need to deal with mean!
