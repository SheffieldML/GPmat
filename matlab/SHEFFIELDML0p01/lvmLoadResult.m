function [model, lbls] = lvmLoadResult(modelType, dataSet, number)

% LVMLOADRESULT Load a previously saved result.
%
%	Description:
%
%	[MODEL, LBLS] = LVMLOADRESULT(DATASET, NUMBER) loads a previously
%	saved LVM result.
%	 Returns:
%	  MODEL - the saved model.
%	  LBLS - labels of the data set (for visualisation purposes).
%	 Arguments:
%	  DATASET - the name of the data set to load.
%	  NUMBER - the number of the LVM run to load.
%	
%
%	See also
%	LVMLOADDATA


%	Copyright (c) 2003, 2004, 2005, 2006, 2008, 2009 Neil D. Lawrence


  fhandle = [modelType 'LoadResult'];
  if exist(fhandle) == 2
    fhandle = str2func(fhandle);
    [model, lbls] = fhandle(dataSet, number);
  else
    [Y, lbls] = lvmLoadData(dataSet);
    
    dataSet(1) = upper(dataSet(1));

    if ~isempty(modelType)
      modelType(1) = upper(modelType(1));
    end

    load(['dem' dataSet modelType num2str(number)])
  end
end
