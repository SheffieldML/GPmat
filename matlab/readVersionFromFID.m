function val = readVersionFromFID(FID)
  
% READDOUBLEFROMFID Read version number from an FID.
% FORMAT
% DESC reads version number from a stream.
% ARG FID : stream to read from.
% RETURN val : value of variable in file.
%
% SEEALSO : writeVersionToFID, readDoubleFromFID, readStringFromFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% NDLUTILS
  
val = str2num(readStringFromFID(FID, 'version'));
