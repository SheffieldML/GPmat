function params = fgplvmExtractParam(model)

% FGPLVMEXTRACTPARAM Extract a parameter vector from a GP-LVM model.
%
% params = fgplvmExtractParam(model)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmExtractParam.m version 1.2



params = gpExtractParam(model);
if isfield(model, 'back')
  params = [modelExtractParam(model.back) params];
else
  params = [model.X(:)' params];
end
if isfield(model, 'dynamics') 
  params = [params modelExtractParam(model.dynamics)];
end