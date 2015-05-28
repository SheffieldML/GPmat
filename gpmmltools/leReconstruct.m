function model = leReconstruct(leInfo, y)

% LERECONSTRUCT Reconstruct an LE form component parts.
% FORMAT
% DESC takes component parts of an LE model and reconstructs the
% LE model. The component parts are normally retrieved from a
% saved file.
% ARG leInfo : the active set and other information stored in a structure.
% ARG y : the output target training data for the LE.
% RETURN model : an LE model structure that combines the component
% parts.
% 
% SEEALSO : leDeconstruct, leCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  model = leInfo;
  model.Y = y;
    
end
