function writeStringToFID(FID, name, val)
  
% WRITESTRINGTOFID Writes a string to an FID.
% FORMAT
% DESC writes an string from a stream.
% ARG FID : stream to write from.
% ARG name : name of string.
% ARG val : value of variable to place in file.
%
% SEEALSO : readStringFromFID, writeIntToFID, writeBoolToFID, writeDoubleToFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% NDLUTIL
  
fprintf(FID, [name '=' val '\n']);
