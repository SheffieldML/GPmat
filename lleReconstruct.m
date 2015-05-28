function model = lleReconstruct(lleInfo, y)

% LLERECONSTRUCT Reconstruct an LLE form component parts.
% FORMAT
% DESC takes component parts of an LLE model and reconstructs the
% LLE model. The component parts are normally retrieved from a
% saved file.
% ARG lleInfo : the active set and other information stored in a structure.
% ARG y : the output target training data for the LLE.
% RETURN model : an LLE model structure that combines the component
% parts.
% 
% SEEALSO : lleDeconstruct, lleCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  model = lleInfo;
  model.Y = y;
    
end
