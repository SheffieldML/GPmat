function model = mapmodelReadFromFID(FID, varargin)

% MAPMODELREADFROMFID Load from a FID produced by C++ code.
% FORMAT
% DESC loads in from a file stream the data format produced by C++ code.
% ARG FID : the file ID from where the data is loaded.
% RETURN model : the model loaded in from the file.
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : modelReadFromFID

% MLTOOLS

modelType = readStringFromFID(FID, 'type');
feval = str2func([modelType 'ReadFromFID']);
model = feval(FID, varargin{:});
