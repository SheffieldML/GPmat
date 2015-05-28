function vec = readBinaryDoubles(fileName, format)

% READBINARYDOUBLES Read information from a binary file in as doubles.
%
%	Description:
%
%	VEC = READBINARYDOUBLES(FILENAME, FORMAT) reads in information from
%	a binary file as a vector of doubles.
%	 Returns:
%	  VEC - vector of values in the file.
%	 Arguments:
%	  FILENAME - the name of the file.
%	  FORMAT - the file format for 'fopen', defaults to 'ieee-le'.
%	
%
%	See also
%	FOPEN, FREAD, FCLOSE


%	Copyright (c) 2009 Neil D. Lawrence

  
  if nargin < 2
    format = 'ieee-le';
  end
  fid = fopen(fileName, 'r', format);
  vec = fread(fid, inf, 'double')';
  fclose(fid);
end