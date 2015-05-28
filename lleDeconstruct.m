function lleInfo = lleDeconstruct(model)

% LLEDECONSTRUCT break LLE in pieces for saving.
% FORMAT
% DESC takes an LLE model structure and breaks it into component
% parts for saving. 
% ARG model : the model that needs to be saved.
% RETURN lleInfo : a structure containing the other information
% from the LLE: what the sparse approximation is, what the inducing
% variables are.
%
% SEEALSO : lleReconstruct
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

lleInfo = model;
removeFields = {'Y'};

for i = 1:length(removeFields)
  if isfield(lleInfo, removeFields{i})
    lleInfo = rmfield(lleInfo, removeFields{i});
  end
end

