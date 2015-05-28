function [model,labels] = fgplvmReadFromFile(fileName)

% FGPLVMREADFROMFILE Load a file produced by the C++ implementation.
% FORMAT
% DESC loads a GP-LVM model from a file produced by the C++
% implementation of the GP-LVM.
% ARG fileName : the file name written by the C++ software.
% RETURN model : a MATLAB GP-LVM model structure containing the
% model from the file.
% RETURN labels : any available labels from the GP-LVM file.
%
% SEEALSO : fgplvmReadFromFID, fgplvmCreate
%
% COPYRIGHT : Neil D. Lawrence, 2005

% FGPLVM

FID = fopen(fileName);
if FID==-1
  error(['Cannot find file ' fileName])
end
[model, labels] = fgplvmReadFromFID(FID);
fclose(FID);
