function model = isomapReconstruct(isomapInfo, y)

% ISOMAPRECONSTRUCT Reconstruct an isomap form component parts.
% FORMAT
% DESC takes component parts of an isomap model and reconstructs the
% isomap model. The component parts are normally retrieved from a
% saved file.
% ARG isomapInfo : the active set and other information stored in a structure.
% ARG y : the output target training data for the isomap.
% RETURN model : an isomap model structure that combines the component
% parts.
% 
% SEEALSO : isomapDeconstruct, isomapCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  model = isomapInfo;
  model.Y = y;
    
end
