function [model, lbls] = modelReadFromFile(fileName, varargin)

% MODELREADFROMFILE Read model from a file FID produced by the C++ implementation.
%
%	Description:
%
%	MODEL = MODELREADFROMFILE(FILENAME) loads in from a file a model
%	produced by C++ code. C++ GP implementation.
%	 Returns:
%	  MODEL - the model loaded in from the file.
%	 Arguments:
%	  FILENAME - the file ID from where the data is loaded.
%	
%
%	See also
%	MODELREADFROMFILE


%	Copyright (c) 2008 Neil D. Lawrence



FID = fopen(fileName);
if FID==-1
  error(['Cannot find file ' fileName])
end
model = modelReadFromFID(FID, varargin{:});
fclose(FID);
