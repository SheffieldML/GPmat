function model = gpReadFromFile(fileName, varargin)

% GPREADFROMFILE Load a file produced by the C++ implementation.
%
%	Description:
%
%	MODEL = GPREADFROMFILE(FILENAME) loads a GP model from a file
%	produced by the C++ GP implementation.
%	 Returns:
%	  MODEL - a MATLAB GP model structure containing the model from the
%	   file.
%	 Arguments:
%	  FILENAME - the file name written by the C++ software.
%	
%
%	See also
%	GPREADFROMFID, GPCREATE


%	Copyright (c) 2005 Neil D. Lawrence


FID = fopen(fileName);
if FID==-1
  error(['Cannot find file ' fileName])
end
model = gpReadFromFID(FID, varargin{:});
fclose(FID);