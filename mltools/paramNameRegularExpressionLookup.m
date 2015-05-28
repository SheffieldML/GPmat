function ind = paramNameRegularExpressionLookup(model, pattern, isKernel)

% PARAMNAMEREGULAREXPRESSIONLOOKUP Returns the indices of the parameter containing the given regular expression.
% FORMAT
% DESC returns the index of the parameters which contain the given
% regular expression. If no matches
% are found then an empty vector is returned.
% ARG model : the model for which parameters are reverse looked up.
% ARG pattern : the regular expression that should match the names.
% RETURN ind : the indices of those parameters in the model.
%
% FORMAT
% DESC returns the index of the parameters which contain the given
% regular expression for a kernel structure. If no matches
% are found then an empty vector is returned.
% ARG model : the model for which parameters are reverse looked up.
% ARG pattern : the regular expression that should match the names.
% ARG isKernel : a logical value that specifies if a model corresponds to a
% kernel structure.
% RETURN ind : the indices of those parameters in the model.
%
% COPYRIGHT : Neil D. Lawrence, 2008
% 
% MODIFICATIONS : Mauricio A. Alvarez, 2009
%
% SEEALSO : cmpndTieParameters, paramNameReverseLookup

% MLTOOLS

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
