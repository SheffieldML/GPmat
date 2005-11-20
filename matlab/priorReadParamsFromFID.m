function prior = priorReadParamsFromFID(prior, FID)

% PRIORREADPARAMSFROMFID Read prior params from C++ written FID.

% PRIOR

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'numParams'))
  error('Incorrect file format.')
end
numParams = str2num(tokens{2});

params = fscanf(FID, '%f\n', numParams);
fhandle = str2func([prior.type 'PriorExpandParam']);
prior = fhandle(prior, params);

