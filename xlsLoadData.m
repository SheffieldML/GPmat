function [data, names] = xlsLoadData(fileName, sheetName)

% XLSLOADDATA Wrapper function for xlsread to get files from the datasets directory.

% GPMAT

baseDir = datasetsDirectory;
[data, names] = xlsread([baseDir fileName], sheetName);
