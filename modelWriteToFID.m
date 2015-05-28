function modelWriteToFID(FID, model)

% MODELWRITETOFID Write to a stream a given model.
% FORMAT
% DESC loads in from a file stream the data format produced by C++ code.
% ARG FID : the file ID from where the data is loaded.
% ARG model : the model loaded in from the file.
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : modelReadFromFID

% MLTOOLS

writeVersionToFID(FID, 0.2);
modelType = readStringFromFID(FID, 'baseType');
feval = str2func([modelType 'ReadFromFID']);
model = feval(FID, varargin{:});
