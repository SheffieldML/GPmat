function model = ppcaReconstruct(ppcaInfo, y)

% PPCARECONSTRUCT Reconstruct an PPCA form component parts.
% FORMAT
% DESC takes component parts of an PPCA model and reconstructs the
% PPCA model. The component parts are normally retrieved from a
% saved file.
% ARG ppcaInfo : the active set and other information stored in a structure.
% ARG y : the output target training data for the PPCA.
% RETURN model : an PPCA model structure that combines the component
% parts.
% 
% SEEALSO : ppcaDeconstruct, ppcaCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  model = ppcaInfo;
  model.y = y;
    
end
