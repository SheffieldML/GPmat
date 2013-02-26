function mvuInfo = mvuDeconstruct(model)

% MVUDECONSTRUCT break MVU in pieces for saving.
%
%	Description:
%
%	MVUINFO = MVUDECONSTRUCT(MODEL) takes an MVU model structure and
%	breaks it into component parts for saving.
%	 Returns:
%	  MVUINFO - a structure containing the other information from the
%	   MVU: what the sparse approximation is, what the inducing variables
%	   are.
%	 Arguments:
%	  MODEL - the model that needs to be saved.
%	
%
%	See also
%	MVURECONSTRUCT


%	Copyright (c) 2009 Neil D. Lawrence


mvuInfo = model;
removeFields = {'Y'};

for i = 1:length(removeFields)
  if isfield(mvuInfo, removeFields{i})
    mvuInfo = rmfield(mvuInfo, removeFields{i});
  end
end

