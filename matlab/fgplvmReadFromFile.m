function [model,labels] = fgplvmReadFromFile(fileName)

% FGPLVMREADFROMFILE Load a file produced by the c++ implementation.

% FGPLVM

FID = fopen(fileName);
if FID==-1
  error(['Cannot find file ' fileName])
end
[model, labels] = fgplvmReadFromFID(FID);
fclose(FID);