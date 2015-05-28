function [model, lbls] = fgplvmLoadResult(dataSet, number)

% FGPLVMLOADRESULT Load a previously saved result.
% FORMAT
% DESC loads a previously saved FGPLVM result.
% ARG dataSet : the name of the data set to load.
% ARG number : the number of the FGPLVM run to load.
% RETURN model : the saved model.
% RETURN lbls : labels of the data set (for visualisation purposes).
%
% SEEALSO : fgplvmLoadResult
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2005, 2006, 2008
  
% FGPLVM

[Y, lbls] = lvmLoadData(dataSet);

dataSet(1) = upper(dataSet(1));
load(['dem' dataSet 'Fgplvm' num2str(number)])
model = fgplvmReconstruct(kern, noise, fgplvmInfo, X, Y);
