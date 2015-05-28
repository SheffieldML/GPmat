function g = cgpsimLogLikeGradients(model)

% CGPSIMLOGLIKEGRADIENTS Compound GPSIM model's gradients.
% FORMAT
% DESC this is just a wrapper for learning three GP models
% simultaneously.
% ARG model : input model for which gradients of the likelihood will be computed.
% RETURN g : gradients of the log likelihood with repsect to the parameters.
%
% SEEALSO : gpsimLogLikeGradients
%
% COPYRIGHT : Neil D. Lawrence, 2006

% SHEFFIELDML

g = gpsimLogLikeGradients(model.comp{1});
for i = 2:length(model.comp)
  g = g + gpsimLogLikeGradients(model.comp{i});
end
