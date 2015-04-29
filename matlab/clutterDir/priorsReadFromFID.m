function priors = priorsReadFromFID(FID, numPriors)

% PRIORSREADFROMFID Read a prior from a C++ written FID.

lineStr = getline(FID);
tokens = tokenise(lineStr, '=')
for i=1:numPriors
  if(length(tokens)~=2 || ~strcmp(tokens{1}, 'priorIndex'))
    error('Incorrect file format.')
  end
  priorIndex = str2num(tokens{2});
  priors{priorIndex}=distReadFromFID(FID);
end
