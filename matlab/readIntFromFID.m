function val = readIntFromFID(FID, string)
  
% READINTFROMFID Read an integer from an FID.
% FORMAT
% DESC reads an integer from a stream.
% ARG FID : stream to read from.
% ARG name : name of integer.
% RETURN val : value of variable in file.
%
% SEEALSO : writeIntToFID, readDoubleFromFID, readIntFromFID,
% readBoolFromFID, readStringFromFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% NDLUTIL
  
  
val = str2num(readStringFromFID(FID, string));
