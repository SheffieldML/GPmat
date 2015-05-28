function ll = cgpsimLogLikelihood(model)

% CGPSIMLOGLIKELIHOOD Compound GPSIM model's log likelihood.
% FORMAT
% DESC this is just a wrapper for learning three GP models
% simultaneously.
% ARG model : input model for which gradients of the likelihood will be computed.
% RETURN ll : the log likelihood.
%
% SEEALSO : gpsimLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2006

% SHEFFIELDML

ll = 0;
for i = 1:length(model.comp)
  ll = ll + gpsimLogLikelihood(model.comp{i});
end
