function gpWriteToFile(fileName, model)

% GPWRITETOFILE Write a file to be read by the C++ implementation.
% FORMAT
% DESC writes a GP model to a file produced by the C++
% GP implementation.
% ARG fileName : the file name written by the C++ software.
% ARG model : a MATLAB GP model structure containing the
% model for the file.
%
% SEEALSO : gpWriteToFID, gpCreate
%
% COPYRIGHT : Neil D. Lawrence, 2008

% GP

FID = fopen(fileName, 'w');
if FID==-1
  error(['Cannot open file ' fileName])
end
gpWriteToFID(FID, model);
fclose(FID);
