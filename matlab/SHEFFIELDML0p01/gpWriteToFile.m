function gpWriteToFile(fileName, model)

% GPWRITETOFILE Write a file to be read by the C++ implementation.
%
%	Description:
%
%	GPWRITETOFILE(FILENAME, MODEL) writes a GP model to a file produced
%	by the C++ GP implementation.
%	 Arguments:
%	  FILENAME - the file name written by the C++ software.
%	  MODEL - a MATLAB GP model structure containing the model for the
%	   file.
%	
%
%	See also
%	GPWRITETOFID, GPCREATE


%	Copyright (c) 2008 Neil D. Lawrence


FID = fopen(fileName, 'w');
if FID==-1
  error(['Cannot open file ' fileName])
end
gpWriteToFID(FID, model);
fclose(FID);