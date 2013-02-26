function [model, lbls] = modelReadFromFID(FID, varargin)

% MODELREADFROMFID Load from a FID produced by C++ code.
%
%	Description:
%
%	MODEL = MODELREADFROMFID(FID) loads in from a file stream the data
%	format produced by C++ code.
%	 Returns:
%	  MODEL - the model loaded in from the file.
%	 Arguments:
%	  FID - the file ID from where the data is loaded.
%	
%
%	See also
%	MODELREADFROMFILE


%	Copyright (c) 2008 Neil D. Lawrence

  
  version = readVersionFromFID(FID);
  if version < 0.2
    error('Incorrect file version.')
  end
  
  modelType = readStringFromFID(FID, 'baseType');
  feval = str2func([modelType 'ReadFromFID']);
  model = feval(FID, varargin{:});
