function params = fgplvmExtractParam(model)

% FGPLVMEXTRACTPARAM Extract a parameter vector from a GP-LVM model.
% FORMAT
% DESC extracts a parameter vector from a given FGPLVM structure.
% ARG model : the model from which parameters are to be extracted.
% RETURN params : the parameter vector extracted from the model.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% SEEALSO : fgplvmCreate, fgplvmExpandParam, modelExtractParam

% FGPLVM

params = gpExtractParam(model);
if isfield(model, 'back')
  params = [modelExtractParam(model.back) params];
else
  params = [model.X(:)' params];
end
if isfield(model, 'dynamics') 
  params = [params modelExtractParam(model.dynamics)];
end