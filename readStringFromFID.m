function val = readStringFromFID(FID, string)
  
% READSTRINGFROMFID Read an boolean from an FID.
% FORMAT
% DESC reads a string from a stream.
% ARG FID : stream to read from.
% ARG name : name of string.
% RETURN val : value of variable in file.
%
% SEEALSO : writeStringToFID, readIntFromFID, readBoolFromFID, readDoubleFromFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% NDLUTIL
  
lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, string))
  error('Incorrect file format.')
end
val = tokens{2};
