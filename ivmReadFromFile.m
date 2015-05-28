function [model, I] = ivmReadFromFile(fileName)

% IVMREADFROMFILE Load a file produced by the C++ implementation.
% FORMAT
% DESC loads an IVM model from a file produced by the C++
% implementation of the IVM.
% ARG fileName : the file name written by the C++ software.
% RETURN model : a MATLAB IVM model structure containing the
% model from the file.
% RETURN I : the indices of the active points in the original data set.
%
% SEEALSO : ivmReadFromFID
%
% COPYRIGHT : Neil D. Lawrence, 2007

% IVM

FID = fopen(fileName);
if FID==-1
  error(['Cannot find file ' fileName])
end
[model, I] = ivmReadFromFID(FID);
fclose(FID);
