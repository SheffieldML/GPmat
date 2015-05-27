function ind = paramNameReverseLookup(model, name)

% PARAMNAMEREVERSELOOKUP Returns the index of the parameter with the given name.
% FORMAT
% DESC returns the index of the parameter with the given name. If no matches
% are found then an empty vector is returned.
% ARG model : the model for which parameters are reverse looked up.
% ARG name : the name of the parameter that is sought.
% RETURN ind : the indices of those parameters in the model.
%
% COPYRIGHT : Neil D. Lawrence, 2008
% 
% SEEALSO : cmpndTieParameters, paramNameRegularExpressionLookup

% MLTOOLS
  
  ind = [];
  [void, names] = modelExtractParam(model);
  for i = 1:length(names)
    if(strcmp(names{i}, name))
      ind = i;
      return
    end
  end
end
