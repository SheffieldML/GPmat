function prior = priorReadFromFID(FID)

% PRIORREADFROMFID Read a prior from a C++ written FID.

% PRIOR

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'distVersion'))
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
prior = priorCreate(type);

prior = priorReadParamsFromFID(prior, FID);


