function model = cmpndTieParameters(model, paramsList)

% CMPNDTIEPARAMETERS Tie parameters together.

% OPTIMI

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
end
