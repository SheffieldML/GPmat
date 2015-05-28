function isomapInfo = isomapDeconstruct(model)

% ISOMAPDECONSTRUCT break isomap in pieces for saving.
% FORMAT
% DESC takes an isomap model structure and breaks it into component
% parts for saving. 
% ARG model : the model that needs to be saved.
% RETURN isomapInfo : a structure containing the other information
% from the isomap: what the sparse approximation is, what the inducing
% variables are.
%
% SEEALSO : isomapReconstruct
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

isomapInfo = model;
removeFields = {'Y'};

for i = 1:length(removeFields)
  if isfield(isomapInfo, removeFields{i})
    isomapInfo = rmfield(isomapInfo, removeFields{i});
  end
end

