function [params, names] = fgplvmExtractParam(model)

% FGPLVMEXTRACTPARAM Extract a parameter vector from a GP-LVM model.
% FORMAT
% DESC extracts a parameter vector from a given FGPLVM structure.
% ARG model : the model from which parameters are to be extracted.
% RETURN params : the parameter vector extracted from the model.
%
% DESC does the same as above, but also returns parameter names.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
% RETURN names : cell array of parameter names.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% SEEALSO : fgplvmCreate, fgplvmExpandParam, modelExtractParam

% FGPLVM

if nargout > 1
  returnNames = true;
else
  returnNames = false;
end  

if returnNames
  [params, names] = gpExtractParam(model);
else
  params = gpExtractParam(model);
end
if isfield(model, 'back') & ~isempty(model.back)
  
  



  if returnNames
    [backParams, backNames] = modelExtractParam(model.back);
    for i = 1:length(backNames)
      backNames{i} = ['Back constraint, ' backNames{i}];
    end
    names = {backNames{:}, names{:}};
  else
    backParams = modelExtractParam(model.back);
  end
  params = [backParams params];
else
  params = [model.X(:)' params];
  if returnNames
    for i = 1:size(model.X, 1)
      for j = 1:size(model.X, 2)
        Xnames{i, j} = ['X(' num2str(i) ', ' num2str(j) ')'];
      end
    end
    names = {Xnames{:}, names{:}};
  end
end
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
  if returnNames
    [dynParams, dynNames] = modelExtractParam(model.dynamics);
    for i = 1:length(dynNames)
      dynNames{i} = ['Dynamics, ' dynNames{i}];
    end
    names = {names{:}, dynNames{:}};
  else
    dynParams = modelExtractParam(model.dynamics);
  end
  params = [params dynParams];
    
end
