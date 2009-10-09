function pmvuInfo = pmvuDeconstruct(model)

% PMVUDECONSTRUCT break PMVU in pieces for saving.
% FORMAT
% DESC takes an probabilistic maximum variance unfolding model structure and breaks it into component
% parts for saving. 
% ARG model : the model that needs to be saved.
% RETURN pmvuInfo : a structure containing the other information
% from the probabilistic maximum variance unfolding: what the sparse approximation is, what the inducing
% variables are.
%
% SEEALSO : pmvuCreate, pmvuReconstruct
%
% COPYRIGHT : Neil D. Lawrence 2009
 
% MLTOOLS


  pmvuInfo = model;
  removeFields = {'Y'};
  
  for i = 1:length(removeFields)
    if isfield(pmvuInfo, removeFields{i})
      pmvuInfo = rmfield(pmvuInfo, removeFields{i});
    end
  end
end
