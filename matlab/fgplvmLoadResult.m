function [model, lbls] = fgplvmLoadResult(dataSet, number)

% FGPLVMLOADRESULT Load a previously saved result.
%
% [model, lbls] = fgplvmLoadResult(dataSet, number)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmLoadResult.m version 1.1



[Y, lbls] = lvmLoadData(dataSet);

dataSet(1) = upper(dataSet(1));
load(['dem' dataSet num2str(number)])
