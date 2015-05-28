function ind = paramNameReverseLookup(model, name)

% PARAMNAMEREVERSELOOKUP Returns the index of the parameter with the given name.
%
%	Description:
%
%	IND = PARAMNAMEREVERSELOOKUP(MODEL, NAME) returns the index of the
%	parameter with the given name. If no matches are found then an empty
%	vector is returned.
%	 Returns:
%	  IND - the indices of those parameters in the model.
%	 Arguments:
%	  MODEL - the model for which parameters are reverse looked up.
%	  NAME - the name of the parameter that is sought.
%	
%
%	See also
%	CMPNDTIEPARAMETERS, PARAMNAMEREGULAREXPRESSIONLOOKUP


%	Copyright (c) 2008 Neil D. Lawrence

  
  ind = [];
  [void, names] = modelExtractParam(model);
  for i = 1:length(names)
    if(strcmp(names{i}, name))
      ind = i;
      return
    end
  end
end