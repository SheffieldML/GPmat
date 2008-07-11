function writeDoubleToFID(FID, name, val)
  
% WRITEDOUBLETOFID Writes a double to an FID.
% FORMAT
% DESC writes a double to a stream.
% ARG FID : stream to write to.
% ARG name : name of double.
% ARG val : value of variable to place in file.
%
% SEEALSO : readDoubleFromFID, writeStringToFID, writeBoolToFID, writeIntToFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% NDLUTIL
  
writeStringToFID(FID, name, num2str(val));
