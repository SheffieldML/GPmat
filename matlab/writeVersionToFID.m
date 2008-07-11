function writeVersionToFID(FID, val)
  
% WRITEVERSIONTOFID Writes a version to an FID.
% FORMAT
% DESC writes a version from a stream.
% ARG FID : stream to write to.
% ARG val : value of version to place in file.
%
% SEEALSO : readVersionFromFID, writeStringToFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% NDLUTIL
  
writeStringToFID(FID, 'version', num2str(val));
  
