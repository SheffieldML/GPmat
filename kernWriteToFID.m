function kernWriteToFID(kern, FID)

% KERNWRITETOFID Load from an FID written by the C++ implementation.
% FORMAT
% DESC loads in from a file stream the data format produced by 
% C++ implementations.
% ARG kern : the kernel structure to write to the stream.
% ARG FID : the file ID from where the data is loaded.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2008
%
% SEEALSO : modelReadFromFID, kernCreate, kernReadParamsFromFID

% KERN

writeVersionToFID(FID, 0.2);
writeStringToFID(FID, 'baseType', 'kern');  
writeStringToFID(FID, 'type', kern.type);
kernWriteParamsToFID(kern, FID);
