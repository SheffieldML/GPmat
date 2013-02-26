function leInfo = leDeconstruct(model)

% LEDECONSTRUCT break LE in pieces for saving.
%
%	Description:
%
%	LEINFO = LEDECONSTRUCT(MODEL) takes an LE model structure and breaks
%	it into component parts for saving.
%	 Returns:
%	  LEINFO - a structure containing the other information from the LE:
%	   what the sparse approximation is, what the inducing variables are.
%	 Arguments:
%	  MODEL - the model that needs to be saved.
%	
%
%	See also
%	LERECONSTRUCT


%	Copyright (c) 2009 Neil D. Lawrence


  leInfo = model;
  removeFields = {'Y'};
  
  for i = 1:length(removeFields)
    if isfield(leInfo, removeFields{i})
      leInfo = rmfield(leInfo, removeFields{i});
    end
  end
end
