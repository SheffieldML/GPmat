function mvuInfo = mvuDeconstruct(model)

% MVUDECONSTRUCT break MVU in pieces for saving.
% FORMAT
% DESC takes an MVU model structure and breaks it into component
% parts for saving. 
% ARG model : the model that needs to be saved.
% RETURN mvuInfo : a structure containing the other information
% from the MVU: what the sparse approximation is, what the inducing
% variables are.
%
% SEEALSO : mvuReconstruct
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

mvuInfo = model;
removeFields = {'Y'};

for i = 1:length(removeFields)
  if isfield(mvuInfo, removeFields{i})
    mvuInfo = rmfield(mvuInfo, removeFields{i});
  end
end

