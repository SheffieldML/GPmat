function [model, lbls] = fgplvmLoadResult(dataSet, number)

% FGPLVMLOADRESULT Load a previously saved result.

% FGPLVM

[Y, lbls] = lvmLoadData(dataSet);

dataSet(1) = upper(dataSet(1));
load(['dem' dataSet num2str(number)])
