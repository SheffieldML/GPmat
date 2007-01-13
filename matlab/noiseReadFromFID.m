function noise = noiseReadFromFID(FID)

% NOISEREADFROMFID Load from an FID written by the C++ implementation.
% FORMAT
% DESC reads a noise structure from a file stream created by the
% C++ implementation of the noises.
% ARG FID : the file stream ID.
% ARG noise : the noise structure created.
%
% COPYRIGHT : Neil D. Lawrence, 2005
%
% SEEALSO : noiseCreate, noiseReadParamsFromFID

% NOISE

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'noiseVersion'))
  error('Incorrect file format.')
end
if(~strcmp(tokens{2}, '0.1'))
  error('Incorrect file version.')
end

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'type'))
  error('Incorrect file format.')
end
type = tokens{2};
noise = noiseReadParamsFromFID(type, FID);

