function params = fgplvmExtractParam(model)

% FGPLVMEXTRACTPARAM Extract a parameter vector from a GP-LVM model.

% FGPLVM

params = gpExtractParam(model);
if isfield(model, 'back')
  params = [modelExtractParam(model.back) params];
else
  params = [model.X(:)' params];
end

