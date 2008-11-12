function model = modelTieParam(model, paramsList)

% MODELTIEPARAM Tie parameters of a model together.
% FORMAT
% DESC groups of parameters of a model to be seen as one parameter
% during optimisation of the model.
% ARG model : the model for which parameters are being tied
% together.
% ARG paramsList : indices of parameteres to group together. The
% indices are provided in a cell array. Each cell in the array
% contains a vector of indices of parameters that should be
% considered as one parameter. Each group of parameters in each
% cell should obviously be mutually exclusive.
% RETURN model : the model with the parameters grouped together.
%
% SEEALSO : modelExtractParam, modeExpandParam, modelLogLikeGradients
% 
% COPYRIGHT : Neil D. Lawrence, 2003, 2006, 2008

% MLTOOLS

  if ~isfield(model, 'paramGroups')
    if isfield(model, 'nParams')
      model.paramGroups = speye(model.nParams);
    elseif isfield(model, 'numParams')
      model.paramGroups = speye(model.numParams);
    else
      error('Model does not list number of parameters.');
    end
  end
  colToDelete = [];
  for i = 1:length(paramsList)
    
    paramIndices=sort(paramsList{i});
    if any(paramIndices(1) == colToDelete)
      error('Parameter is already being tied')
    end
    for j = 2:length(paramIndices)
      
      model.paramGroups(paramIndices(j), paramIndices(1)) = 1;
      if any(paramIndices(j) == colToDelete)
        error('Parameter has already been tied')
      end
      colToDelete = [colToDelete paramIndices(j)];
    end
  end
  
  model.paramGroups(:, colToDelete) = [];
  if isfield(model, 'nParams')
    % Update to the new number of parameters.
    model.nParams = size(model.paramGroups, 2);
  elseif isfield(model, 'numParams')
    model.numParams = size(model.paramGroups, 2);
  end
end