function kern = cmpndTieParameters(kern, paramsList)

% CMPNDTIEPARAMETERS Tie parameters together.

% IVM


colToDelete = [];
for i = 1:length(paramsList)

  paramIndices=sort(paramsList{i});
  if any(paramIndices(1) == colToDelete)
    error('Parameter is already being tied')
  end
  for j = 2:length(paramIndices)

    kern.paramGroups(paramIndices(j), paramIndices(1)) = 1;
    if any(paramIndices(j) == colToDelete)
      error('Parameter has already been tied')
    end
    colToDelete = [colToDelete paramIndices(j)];
  end
end

kern.paramGroups(:, colToDelete) = [];
