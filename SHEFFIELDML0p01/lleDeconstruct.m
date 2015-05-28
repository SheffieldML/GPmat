function lleInfo = lleDeconstruct(model)

% LLEDECONSTRUCT break LLE in pieces for saving.
%
%	Description:
%
%	LLEINFO = LLEDECONSTRUCT(MODEL) takes an LLE model structure and
%	breaks it into component parts for saving.
%	 Returns:
%	  LLEINFO - a structure containing the other information from the
%	   LLE: what the sparse approximation is, what the inducing variables
%	   are.
%	 Arguments:
%	  MODEL - the model that needs to be saved.
%	
%
%	See also
%	LLERECONSTRUCT


%	Copyright (c) 2009 Neil D. Lawrence


lleInfo = model;
removeFields = {'Y'};

for i = 1:length(removeFields)
  if isfield(lleInfo, removeFields{i})
    lleInfo = rmfield(lleInfo, removeFields{i});
  end
end

