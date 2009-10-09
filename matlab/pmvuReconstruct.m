function model = pmvuReconstruct(pmvuInfo, y)

% PMVURECONSTRUCT Reconstruct an probabilistic maximum variance unfolding from component parts.
% FORMAT
% DESC takes component parts of an probabilistic maximum variance unfolding model and reconstructs the
% probabilistic maximum variance unfolding model. The component parts are normally retrieved from a
% saved file.
% ARG pmvuInfo : the active set and other information stored in a structure.
% ARG y : the output target training data for the probabilistic maximum variance unfolding.
% RETURN model : an probabilistic maximum variance unfolding model structure that combines the component
% parts.
%
% SEEALSO : pmvuCreate, pmvuReconstruct
%
% COPYRIGHT : Neil D. Lawrence 2009
 
% MLTOOLS

  model = pmvuInfo;
  model.Y = y;
  
end