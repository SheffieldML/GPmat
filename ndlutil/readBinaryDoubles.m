function vec = readBinaryDoubles(fileName, format)

% READBINARYDOUBLES Read information from a binary file in as doubles.
% FORMAT
% DESC reads in information from a binary file as a vector of doubles.
% ARG fileName : the name of the file.
% ARG format : the file format for 'fopen', defaults to 'ieee-le'.
% RETURN vec : vector of values in the file. 
%
% COPYRIGHT : Neil D. Lawrence, 2009
%
% SEEALSO : fopen, fread, fclose

% NDLUTIL
  
  if nargin < 2
    format = 'ieee-le';
  end
  fid = fopen(fileName, 'r', format);
  vec = fread(fid, inf, 'double')';
  fclose(fid);
end
