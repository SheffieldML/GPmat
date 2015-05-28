function [param, names] = cgpdisimExtractParam(model)

% CGPDISIMEXTRACTPARAM Extract parameters from compound GPDISIM model.
% FORMAT
% DESC this is just a wrapper for learning three GP models
% simultaneously.
% ARG model : input model from which parameters will be extracted.
% RETURN param : parameters extracted from the model.
%
% SEEALSO : gpdisimExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007

% SHEFFIELDML

if nargout == 1,
  param = gpdisimExtractParam(model.comp{1});
else
  [param, names] = gpdisimExtractParam(model.comp{1});
end
