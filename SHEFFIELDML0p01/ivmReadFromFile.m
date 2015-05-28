function [model, I] = ivmReadFromFile(fileName)

% IVMREADFROMFILE Load a file produced by the C++ implementation.
%
%	Description:
%
%	[MODEL, I] = IVMREADFROMFILE(FILENAME) loads an IVM model from a
%	file produced by the C++ implementation of the IVM.
%	 Returns:
%	  MODEL - a MATLAB IVM model structure containing the model from the
%	   file.
%	  I - the indices of the active points in the original data set.
%	 Arguments:
%	  FILENAME - the file name written by the C++ software.
%	
%
%	See also
%	IVMREADFROMFID


%	Copyright (c) 2007 Neil D. Lawrence


FID = fopen(fileName);
if FID==-1
  error(['Cannot find file ' fileName])
end
[model, I] = ivmReadFromFID(FID);
fclose(FID);