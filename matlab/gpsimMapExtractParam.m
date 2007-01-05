function [param, names] = gpsimMapExtractParam(model)

% GPSIMMAPEXTRACTPARAM Extract the parameters of a GPSIMMAP model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a Gaussian process for single input motif
% modelling.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% SEEALSO : gpsimMapCreate, gpsimMapExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GPSIM

if nargout>1
  [param, names] = kernExtractParam(model.kern);
else
  param = kernExtractParam(model.kern);
end

if isfield(model, 'fix')
  for i = 1:length(model.fix)
    param(model.fix(i).index) = model.fix(i).value;
  end
end
param = real(param);

% Check if there is a mean function.
if isfield(model, 'meanFunction') & ~isempty(model.meanFunction)
  if nargout>1
    [meanFuncParams, meanFuncNames] = modelExtractParam(model.meanFunction);
  else
    meanFuncParams = modelExtractParam(model.meanFunction);
  end
else
  meanFuncParams =[];
  meanFuncNames = [];
end

param = [param meanFuncParams];
if nargout > 1
  names = {names{:} meanFuncNames{:}};
end