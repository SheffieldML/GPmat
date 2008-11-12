function ind = paramNameRegularExpressionLookup(model, pattern)

% PARAMNAMEREGULAREXPRESSIONLOOKUP Returns the indices of the parameter containing the given regular expression.
% FORMAT
% DESC returns the index of the parameters which contain the given
% regular expression. If no matches
% are found then an empty vector is returned.
% ARG model : the model for which parameters are reverse looked up.
% ARG pattern : the regular expression that should match the names.
% RETURN ind : the indices of those parameters in the model.
%
% COPYRIGHT : Neil D. Lawrence, 2008
% 
% SEEALSO : cmpndTieParameters, paramNameReverseLookup

% MLTOOLS
  
  ind = [];
  [void, names] = modelExtractParam(model);
  for i = 1:length(names)
    if(regexp(names{i}, pattern))
      ind = [ind i];
    end
  end
end