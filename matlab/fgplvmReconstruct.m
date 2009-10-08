function model = fgplvmReconstruct(kern, noise, fgplvmInfo, X, y)

% FGPLVMRECONSTRUCT Reconstruct an FGPLVM from component parts.
% FORMAT
% DESC takes component parts of an FGPLVM model and reconstructs the
% FGPLVM model. The component parts are normally retrieved from a
% saved file.
% ARG kern : a kernel structure for the FGPLVM.
% ARG noise : a noise structure for the FGPLVM.
% ARG fgplvmInfo : the active set and the inactive set of the FGPLVM as
% well as the site parameters, stored in a structure.
% ARG X : the input training data for the FGPLVM.
% ARG y : the output target training data for the FGPLVM.
% RETURN model : an FGPLVM model structure that combines the component
% parts.
% 
% SEEALSO : fgplvmDeconstruct, fgplvmCreate, gpReconstruct
%
% COPYRIGHT : Neil D. Lawrence, 2009

% FGPLVM

  model = gpReconstruct(kern, noise, fgplvmInfo, X, y);
  model.type = 'fgplvm';
  if isfield(model, 'back') && ~isempty(model.back)
    switch model.back.type
     case 'kbr'
      model.back.X = model.y;
     otherwise
      
    end
  end
  params = fgplvmExtractParam(model);
  model = fgplvmExpandParam(model, params);
end
