function val = readBoolFromFID(FID, string)
  
% READBOOLFROMFID Read a boolean from an FID.
% FORMAT
% DESC reads a boolean from a stream.
% ARG FID : stream to read from.
% ARG name : name of boolean.
% RETURN val : value of variable in file.
%
% SEEALSO : writeBoolToFID, readIntFromFID, readDoubleFromFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% NDLUTIL

  
val = str2num(readStringFromFID(FID, string));
