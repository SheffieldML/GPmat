function ll = cgpdisimLogLikelihood(model)

% CGPDISIMLOGLIKELIHOOD Compound GPDISIM model's log likelihood.
% FORMAT
% DESC this is just a wrapper for learning three GP models
% simultaneously.
% ARG model : input model for which gradients of the likelihood will be computed.
% RETURN ll : the log likelihood.
%
% SEEALSO : gpdisimLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007

% SHEFFIELDML

ll = 0;
for i = 1:length(model.comp)
  ll = ll + gpdisimLogLikelihood(model.comp{i});
end
