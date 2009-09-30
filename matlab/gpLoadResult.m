function model = gpLoadResult(dataSet, number)

% GPLOADRESULT Load a previously saved result.
% FORMAT
% DESC loads a previously saved GP result.
% ARG dataSet : the name of the data set to load.
% ARG number : the number of the GP run to load.
% RETURN model : the saved model.
%
% SEEALSO : gpLoadResult
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2005, 2006, 2008
  
% GP

  [X, y] = mapLoadData(dataSet);

  dataSet(1) = upper(dataSet(1));
  load(['dem' dataSet 'Gp' num2str(number)])
  model = gpReconstruct(kern, noise, gpInfo, X, y);
end