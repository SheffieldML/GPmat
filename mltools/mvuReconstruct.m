function model = mvuReconstruct(mvuInfo, y)

% MVURECONSTRUCT Reconstruct an MVU form component parts.
% FORMAT
% DESC takes component parts of an MVU model and reconstructs the
% MVU model. The component parts are normally retrieved from a
% saved file.
% ARG mvuInfo : the active set and other information stored in a structure.
% ARG y : the output target training data for the MVU.
% RETURN model : an MVU model structure that combines the component
% parts.
% 
% SEEALSO : mvuDeconstruct, mvuCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  model = mvuInfo;
  model.Y = y;
    
end
