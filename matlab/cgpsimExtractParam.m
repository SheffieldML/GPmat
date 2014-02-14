function [param, names] = cgpsimExtractParam(model)

% CGPSIMEXTRACTPARAM Extract parameters from compound GPSIM model.
% FORMAT
% DESC this is just a wrapper for learning three GP models
% simultaneously.
% ARG model : input model from which parameters will be extracted.
% RETURN param : parameters extracted from the model.
%
% SEEALSO : gpsimExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% SHEFFIELDML

if nargout == 1,
  param = gpsimExtractParam(model.comp{1});
else
  [param, names] = gpsimExtractParam(model.comp{1});
end
