function kern = kernReadParamsFromFID(kern, FID)

% KERNREADPARAMSFROMFID Read the kernel parameters from C++ file FID.

% KERN

if strcmp(kern.type, 'cmpnd')
  kern = cmpndKernReadParamsFromFID(kern, FID);
else
  lineStr = getline(FID);
  tokens = tokenise(lineStr, '=');
  if(length(tokens)~=2 | ~strcmp(tokens{1}, 'numParams'))
    error('Incorrect file format.')
  end
  numParams = str2num(tokens{2});
  
  lineStr = getline(FID);
  tokens = tokenise(lineStr, '=');
  if(length(tokens)~=2 | ~strcmp(tokens{1}, 'numFeatures'))
    error('Incorrect file format.')
  end
  numFeatures = str2num(tokens{2});

  if strcmp(kern.type, 'poly') | strcmp(kern.type, 'polyard')
    lineStr = getline(FID);
    tokens = tokenise(lineStr, '=');
    if(length(tokens)~=2 | ~strcmp(tokens{1}, 'degree'))
      error('Incorrect file format.')
    end
    kern.degree = str2num(tokens{2});
  end
  
  params = fscanf(FID, '%f\n', numParams);
  kern.inputDimension=numFeatures;
  fhandle = str2func([kern.type 'KernExpandParam']);
  kern = fhandle(kern, params);
  
  lineStr = getline(FID);
  tokens = tokenise(lineStr, '=');
  if(length(tokens)~=2 | ~strcmp(tokens{1}, 'numPriors'))
    error('Incorrect file format.')
  end
  numPriors = str2num(tokens{2});
  for j=1:numPriors
    lineStr = getline(FID);
    tokens = tokenise(lineStr, '=');
    if(length(tokens)~=2 | ~strcmp(tokens{1}, 'priorIndex'))
      error('Incorrect file format.')
    end
    prior = priorReadFromFID(FID);
    prior.index = str2num(tokens{2})+1;
    kern.priors(j) = prior;
  end
end

