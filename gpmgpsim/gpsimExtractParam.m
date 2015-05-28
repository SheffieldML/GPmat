function [param, names] = gpsimExtractParam(model)

% GPSIMEXTRACTPARAM Extract the parameters of a GPSIM model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a Gaussian process for single input motif
% modelling.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% SEEALSO : gpsimCreate, gpsimExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% SHEFFIELDML

if nargout>1
  [param, names] = kernExtractParam(model.kern);
  for i=1:length(model.mu);
    names{end+1}=['Basal transcription ' num2str(i)];
  end
else
  param = kernExtractParam(model.kern);
end
fhandle = str2func([model.bTransform 'Transform']);
param = [param fhandle(model.B, 'xtoa')];

if isfield(model, 'fix')
  for i = 1:length(model.fix)
    param(model.fix(i).index) = model.fix(i).value;
  end
end
param = real(param);
