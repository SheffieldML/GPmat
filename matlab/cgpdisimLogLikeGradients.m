function g = cgpdisimLogLikeGradients(model)

% CGPDISIMLOGLIKEGRADIENTS Compound GPDISIM model's gradients.
% FORMAT
% DESC this is just a wrapper for learning three GP models
% simultaneously.
% ARG model : input model for which gradients of the likelihood will be computed.
% RETURN g : gradients of the log likelihood with repsect to the parameters.
%
% SEEALSO : gpdisimLogLikeGradients
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007

% SHEFFIELDML

g = gpdisimLogLikeGradients(model.comp{1});
for i = 2:length(model.comp)
  g = g + gpdisimLogLikeGradients(model.comp{i});
end
