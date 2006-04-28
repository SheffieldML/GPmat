function [model,labels] = fgplvmReadFromFile(fileName)

% FGPLVMREADFROMFILE Load a file produced by the c++ implementation.
%
% [model,labels] = fgplvmReadFromFile(fileName)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmReadFromFile.m version 1.1



FID = fopen(fileName);
if FID==-1
  error(['Cannot find file ' fileName])
end
[model, labels] = fgplvmReadFromFID(FID);
fclose(FID);