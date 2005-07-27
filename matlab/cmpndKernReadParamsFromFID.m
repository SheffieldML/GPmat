function kern = cmpndKernReadParamsFromFID(kern, FID)

% CMPNDKERNREADPARAMSFROMFID Read a compound kernel from a C++ file.

% KERN

lineStr=getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'numKerns'))
  error('Incorrect file format.')
end
numKerns=str2num(tokens{2});

lineStr=getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'numFeatures'))
  error('Incorrect file format.')
end
kern.inputDimension=str2num(tokens{2});

for i=1:numKerns
  kern.comp{i} = kernReadFromFID(FID);
end

for i = 1:length(kern.comp)
  kern.nParams = kern.nParams + kern.comp{i}.nParams;
  kern.comp{i}.index = [];
end
kern.paramGroups = speye(kern.nParams);

% Summarise the total white variance in the field whiteVariance.
kern.whiteVariance = 0;
for i = 1:length(kern.comp)
  if strcmp(kern.comp{i}.type, 'white')
    kern.whiteVariance = kern.whiteVariance + kern.comp{i}.variance;
  else
    if(isfield(kern.comp{i}, 'whiteVariance'))
      kern.whiteVariance = kern.whiteVariance + ...
          kern.comp{i}.whiteVariance;
    end
  end
end

