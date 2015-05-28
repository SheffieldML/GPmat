function model = gpReconstruct(kern, noise, gpInfo, X, y)

% GPRECONSTRUCT Reconstruct an GP form component parts.
% FORMAT
% DESC takes component parts of an GP model and reconstructs the
% GP model. The component parts are normally retrieved from a
% saved file.
% ARG kern : a kernel structure for the GP.
% ARG noise : a noise structure for the GP (currently ignored).
% ARG gpInfo : the active set and other information stored in a structure.
% ARG X : the input training data for the GP.
% ARG y : the output target training data for the GP.
% RETURN model : an GP model structure that combines the component
% parts.
% 
% SEEALSO : gpDeconstruct, gpCreate
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2009

% GP

  model = gpInfo;
  model.X = X;
  model.y = y;
  model.kern = kern;
  if ~isempty(noise)
    model.noise = noise;
  end
  model.m = gpComputeM(model);
  
  if isfield(model, 'computeS') && model.computeS 
    model.S = model.m*model.m';
  end
  params = gpExtractParam(model);
  model = gpExpandParam(model, params);
  
end
