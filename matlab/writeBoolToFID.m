function writeBoolToFID(FID, name, val)
  
% WRITEBOOLTOFID Writes a boolean to an FID.
% FORMAT
% DESC writes a boolean to a stream.
% ARG FID : stream to write to.
% ARG name : name of boolean.
% ARG val : value of variable to place in file.
%
% SEEALSO : readBoolFromFID, writeStringToFID, writeBoolToFID, writeIntToFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% NDLUTIL
  
if val
  writeStringToFID(FID, name, '1');
else
  writeStringToFID(FID, name, '0');
end
  
