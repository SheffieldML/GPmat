function writeIntToFID(FID, name, val)
  
% WRITEINTTOFID Writes an integer to an FID.
% FORMAT
% DESC writes an integer to a stream.
% ARG FID : stream to write to.
% ARG name : name of int.
% ARG val : value of variable to place in file.
%
% SEEALSO : readIntFromFID, writeStringToFID, writeBoolToFID, writeDoubleToFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% NDLUTIL
  
writeStringToFID(FID, name, num2str(val));
