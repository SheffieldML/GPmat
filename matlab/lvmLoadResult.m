function [model, lbls] = lvmLoadResult(modelType, dataSet, number)

% LVMLOADRESULT Load a previously saved result.
% FORMAT
% DESC loads a previously saved LVM result.
% ARG dataSet : the name of the data set to load.
% ARG number : the number of the LVM run to load.
% RETURN model : the saved model.
% RETURN lbls : labels of the data set (for visualisation purposes).
%
% SEEALSO : lvmLoadData
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2005, 2006, 2008
  
% MLTOOLS



[Y, lbls] = lvmLoadData(dataSet);

dataSet(1) = upper(dataSet(1));
if ~isempty(modelType)
    modelType(1) = upper(modelType(1));
end

load(['dem' dataSet modelType num2str(number)])
