function [model, lbls] = modelReadFromFID(FID, varargin)

% MODELREADFROMFID Load from a FID produced by C++ code.
% FORMAT
% DESC loads in from a file stream the data format produced by C++ code.
% ARG FID : the file ID from where the data is loaded.
% RETURN model : the model loaded in from the file.
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : modelReadFromFile

% MLTOOLS
  
  version = readVersionFromFID(FID);
  if version < 0.2
    error('Incorrect file version.')
  end
  
  modelType = readStringFromFID(FID, 'baseType');
  feval = str2func([modelType 'ReadFromFID']);
  model = feval(FID, varargin{:});
