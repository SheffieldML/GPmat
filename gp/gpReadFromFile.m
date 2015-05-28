function model = gpReadFromFile(fileName, varargin)

% GPREADFROMFILE Load a file produced by the C++ implementation.
% FORMAT
% DESC loads a GP model from a file produced by the C++
% GP implementation.
% ARG fileName : the file name written by the C++ software.
% RETURN model : a MATLAB GP model structure containing the
% model from the file.
%
% SEEALSO : gpReadFromFID, gpCreate
%
% COPYRIGHT : Neil D. Lawrence, 2005

% GP

FID = fopen(fileName);
if FID==-1
  error(['Cannot find file ' fileName])
end
model = gpReadFromFID(FID, varargin{:});
fclose(FID);
