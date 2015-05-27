function g = gpsimMapLogLikeGradients(model, gfdata)

% GPSIMMAPLOGLIKEGRADIENTS Compute the gradients of the log likelihood of a GPSIMMAP model.
% FORMAT
% DESC computes the gradients of the log likelihood of the given
% Gaussian process for use in a single input motif protein network.
% ARG model : the model for which the log likelihood is computed.
% RETURN g : the gradients of the parameters of the model.
%
% FORMAT
% DESC computes the gradients of the log likelihood of the given
% Gaussian process for use in a single input motif protein network.
% ARG model : the model for which the log likelihood is computed.
% ARG gfdata : the functional gradient of the data portion of the log likelihood
% RETURN g : the gradients of the parameters of the model.
% 
% SEEALSO : gpsimMapCreate, gpsimMapLogLikelihood,
% gpsimMapGradient, gpsimMapFunctionalLogLikeGradients
%
% COPYRIGHT : Magnus Rattray and Neil D. Lawrence, 2006
%  
% MODIFIED : Pei Gao, 2008

% SHEFFIELDML

% kernel parameters
dataRespond = false;
if dataRespond & nargin < 2
  [void, gfdata] = gpsimMapFunctionalLogLikeGradients(model);
end
covGrad = model.invK*(model.f*model.f'+model.covf)*model.invK  ...
           - model.invK;
covGrad = covGrad*0.5;
if dataRespond
  covGrad = covGrad - model.invK*model.covf*diag(model.covf* ...
                                                 model.W)*gfdata;
end
gkern = kernGradient(model.kern, model.mapt, covGrad);

% model parameters
gmodel = gpsimMapMarginalLikeliGradient(model);

g = [gkern gmodel];

% Add term arising from data likelihood (pg 125 in Rasmussen & Williams).
if strcmp(model.nonLinearity, 'linear')
  g1 = 0;
else
%  g1 = 0;
  g1 = gpsimMapLikeGradientImplicit(model);
end

g = g + g1;

if isfield(model, 'fix')
  for i = 1:length(model.fix)
    g(model.fix(i).index) = 0;
  end
end
