function leInfo = leDeconstruct(model)

% LEDECONSTRUCT break LE in pieces for saving.
% FORMAT
% DESC takes an LE model structure and breaks it into component
% parts for saving. 
% ARG model : the model that needs to be saved.
% RETURN leInfo : a structure containing the other information
% from the LE: what the sparse approximation is, what the inducing
% variables are.
%
% SEEALSO : leReconstruct
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  leInfo = model;
  removeFields = {'Y'};
  
  for i = 1:length(removeFields)
    if isfield(leInfo, removeFields{i})
      leInfo = rmfield(leInfo, removeFields{i});
    end
  end
end
