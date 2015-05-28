function model = dnetReconstruct(mapping, dnetInfo, y)

% DNETRECONSTRUCT Reconstruct an DNET form component parts.
% FORMAT
% DESC takes component parts of an DNET model and reconstructs the
% DNET model. The component parts are normally retrieved from a
% saved file.
% ARG kern : a kernel structure for the DNET.
% ARG noise : a noise structure for the DNET (currently ignored).
% ARG dnetInfo : the active set and other information stored in a structure.
% ARG X : the input training data for the DNET.
% ARG y : the output target training data for the DNET.
% RETURN model : an DNET model structure that combines the component
% parts.
% 
% SEEALSO : dnetDeconstruct, dnetCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  model = dnetInfo;
  model.mapping = mapping;
  model.y = y;
  
  params = dnetExtractParam(model);
  model = dnetExpandParam(model, params);
  model = dnetEstep(model);
  
end
