function noise = noiseReadFromFID(FID)

% NOISEREADFROMFID Load from an FID written by the C++ implementation.
% FORMAT
% DESC reads a noise structure from a file stream created by the
% C++ implementation of the noises.
% ARG FID : the file stream ID.
% ARG noise : the noise structure created.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2008
%
% SEEALSO : modelReadFromFID, noiseCreate, noiseReadParamsFromFID

% NOISE

type = readStringFromFID(FID, 'type');
noise = noiseReadParamsFromFID(type, FID);

