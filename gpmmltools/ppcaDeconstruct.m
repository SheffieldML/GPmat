function ppcaInfo = ppcaDeconstruct(model)

% PPCADECONSTRUCT break PPCA in pieces for saving.
% FORMAT
% DESC takes an PPCA model structure and breaks it into component
% parts for saving. 
% ARG model : the model that needs to be saved.
% RETURN ppcaInfo : a structure containing the other information
% from the PPCA: what the sparse approximation is, what the inducing
% variables are.
%
% SEEALSO : ppcaReconstruct
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

ppcaInfo = model;
removeFields = {'y'};

for i = 1:length(removeFields)
  if isfield(ppcaInfo, removeFields{i})
    ppcaInfo = rmfield(ppcaInfo, removeFields{i});
  end
end

