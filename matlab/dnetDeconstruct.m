function [mapping, dnetInfo] = dnetDeconstruct(model)

% DNETDECONSTRUCT break DNET in pieces for saving.
% FORMAT
% DESC takes an DNET model structure and breaks it into component
% parts for saving. 
% ARG model : the model that needs to be saved.
% RETURN mapping : the mapping component of the DNET model.
% RETURN dnetInfo : a structure containing the other information
% from the DNET: what the sparse approximation is, what the inducing
% variables are.
%
% SEEALSO : dnetReconstruct
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

mapping = model.mapping;
dnetInfo = model;
removeFields = {'mapping', 'A',  'b', 'Phi', 'y', 'w', 'X'};

for i = 1:length(removeFields)
  if isfield(dnetInfo, removeFields{i})
    dnetInfo = rmfield(dnetInfo, removeFields{i});
  end
end

