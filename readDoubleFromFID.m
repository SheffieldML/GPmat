function val = readDoubleFromFID(FID, string)
  
% READDOUBLEFROMFID Read a double from an FID.
% FORMAT
% DESC reads a double from a stream.
% ARG FID : stream to read from.
% ARG name : name of double.
% RETURN val : value of variable in file.
%
% SEEALSO : writeDoubleToFID, readIntFromFID, readBoolFromFID, readStringFromFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% NDLUTIL
  
val = str2num(readStringFromFID(FID, string));
