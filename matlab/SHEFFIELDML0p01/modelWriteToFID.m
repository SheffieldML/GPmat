function modelWriteToFID(FID, model)

% MODELWRITETOFID Write to a stream a given model.
%
%	Description:
%
%	MODELWRITETOFID(FID, MODEL) loads in from a file stream the data
%	format produced by C++ code.
%	 Arguments:
%	  FID - the file ID from where the data is loaded.
%	  MODEL - the model loaded in from the file.
%	
%
%	See also
%	MODELREADFROMFID


%	Copyright (c) 2008 Neil D. Lawrence


writeVersionToFID(FID, 0.2);
modelType = readStringFromFID(FID, 'baseType');
feval = str2func([modelType 'ReadFromFID']);
model = feval(FID, varargin{:});
