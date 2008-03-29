function kern = kernReadFromFID(FID)

% KERNREADFROMFID Load from an FID written by the C++ implementation.
% FORMAT
% DESC loads in from a file stream the data format produced by 
% C++ implementations.
% ARG FID : the file ID from where the data is loaded.
% RETURN kern : the kernel loaded in from the file.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2008
%
% SEEALSO : modelReadFromFID, kernCreate, kernReadParamsFromFID

% KERN

type = readStringFromFID(FID, 'type');
kern = kernCreate(zeros(1), type);

kern = kernReadParamsFromFID(kern, FID);
