function noiseWriteToFID(noise, FID)

% NOISEWRITETOFID Load from an FID written by the C++ implementation.
% FORMAT
% DESC loads in from a file stream the data format produced by 
% C++ implementations.
% ARG noise : the noise model to write to the stream.
% ARG FID : the file ID from where the data is loaded.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2008
%
% SEEALSO : modelReadFromFID, noiseCreate, noiseReadParamsFromFID

% NOISE

writeVersionToFID(FID, 0.2);
writeStringToFID(FID, 'baseType', 'noise');  
writeStringToFID(FID, 'type', noise.type);
noiseWriteParamsToFID(noise, FID);
