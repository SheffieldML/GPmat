function model = gpLoadResult(dataSet, number, dataLoaderStr)

% GPLOADRESULT Load a previously saved result.
%
%	Description:
%
%	MODEL = GPLOADRESULT(DATASET, NUMBER) loads a previously saved GP
%	result.
%	 Returns:
%	  MODEL - the saved model.
%	 Arguments:
%	  DATASET - the name of the data set to load.
%	  NUMBER - the number of the GP run to load.
%	
%
%	See also
%	GPLOADRESULT


%	Copyright (c) 2003, 2004, 2005, 2006, 2008, 2011 Neil D. Lawrence


  if nargin < 3 
    dataLoaderStr = 'mapLoadData';
  end
  dataLoaderHandler = str2func(dataLoaderStr);
  [X, y] = dataLoaderHandler(dataSet);

  dataSet(1) = upper(dataSet(1));
  load(['dem' dataSet 'Gp' num2str(number)])
  model = gpReconstruct(kern, noise, gpInfo, X, y);
end