function model = cgpdisimExpandParam(model, params)

% CGPDISIMEXPANDPARAM Expand params into model structure.
% FORMAT
% DESC this is just a wrapper for learning three GP models
% simultaneously.
% ARG model : input model from which parameters will be extracted.
% ARG params : parameters to place in the model.
% RETURN model : model filled with the given parameters.
%
% SEEALSO : gpdisimExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007

% SHEFFIELDML

for i = 1:length(model.comp)
  model.comp{i} = gpdisimExpandParam(model.comp{i}, params);
end
