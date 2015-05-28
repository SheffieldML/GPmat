function ind = paramNameRegularExpressionLookup(model, pattern, isKernel)

% PARAMNAMEREGULAREXPRESSIONLOOKUP Returns the indices of the parameter containing the given regular expression.
%
%	Description:
%
%	IND = PARAMNAMEREGULAREXPRESSIONLOOKUP(MODEL, PATTERN) returns the
%	index of the parameters which contain the given regular expression.
%	If no matches are found then an empty vector is returned.
%	 Returns:
%	  IND - the indices of those parameters in the model.
%	 Arguments:
%	  MODEL - the model for which parameters are reverse looked up.
%	  PATTERN - the regular expression that should match the names.
%
%	IND = PARAMNAMEREGULAREXPRESSIONLOOKUP(MODEL, PATTERN, ISKERNEL)
%	returns the index of the parameters which contain the given regular
%	expression for a kernel structure. If no matches are found then an
%	empty vector is returned.
%	 Returns:
%	  IND - the indices of those parameters in the model.
%	 Arguments:
%	  MODEL - the model for which parameters are reverse looked up.
%	  PATTERN - the regular expression that should match the names.
%	  ISKERNEL - a logical value that specifies if a model corresponds
%	   to a kernel structure.
%	
%	
%
%	See also
%	CMPNDTIEPARAMETERS, PARAMNAMEREVERSELOOKUP


%	Copyright (c) 2008 Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009


  ind = [];
  
  if nargin<3
      [void, names] = modelExtractParam(model);
  else
      if isKernel
          [void, names] = kernExtractParam(model);
      else
          [void, names] = modelExtractParam(model);
      end
  end
      
  for i = 1:length(names)
    if(regexp(names{i}, pattern))
      ind = [ind i];
    end
  end
end