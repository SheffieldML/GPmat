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

% GPSIM

dataRespond = false;
if dataRespond & nargin < 2
  [void, gfdata] = gpsimMapFunctionalLogLikeGradients(model);
end
covGrad = model.invK*(model.f*model.f'+model.covf)*model.invK  ...
           - model.invK;
covGrad = covGrad*.5;
if dataRespond
  covGrad = covGrad - model.invK*model.covf*diag(model.covf* ...
                                                 model.W)*gfdata;
end
g = kernGradient(model.kern, model.mapt, covGrad);

% Add term arising from data likelihood (pg 125 in Rasmussen & Williams).
g = g;
