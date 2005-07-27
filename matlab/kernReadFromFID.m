function kern = kernReadFromFID(FID)

% KERNREADFROMFID Load from an FID written by the C++ implementation.

% KERN

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'kernVersion'))
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
kern = kernCreate(zeros(1), type);

kern = kernReadParamsFromFID(kern, FID);
