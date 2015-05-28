function [model, lbls] = modelReadFromFile(fileName, varargin)

% MODELREADFROMFILE Read model from a file FID produced by the C++ implementation.
% FORMAT
% DESC loads in from a file a model produced by C++ code.
% C++ GP implementation.
% ARG fileName : the file ID from where the data is loaded.
% RETURN model : the model loaded in from the file.
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : modelReadFromFile

% MLTOOLS


FID = fopen(fileName);
if FID==-1
  error(['Cannot find file ' fileName])
end
model = modelReadFromFID(FID, varargin{:});
fclose(FID);
