function kern = componentKernReadParamsFromFID(kern, FID, version)

% COMPONENTKERNREADPARAMSFROMFID Read a component based kernel from a C++ file.
% FORMAT
% DESC reads the components fo a kernel from a file written by C++ code.
% ARG kern : the base kernel to add components to.
% ARG FID : the input file stream to add from.
% RETURN kern : the kernel with the components added in.
%
% SEEALSO : modelReadFromFID, kernReadFromFID
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2008
  
% KERN

kern.inputDimension = readIntFromFID(FID, 'inputDim');
numParams = readIntFromFID(FID, 'numParams');
numKerns = readIntFromFID(FID, 'numKerns');

for i=1:numKerns
  if version > 0.11
    kern.comp{i} = modelReadFromFID(FID);
  else
    kern.comp{i} = kernReadFromFID(FID, version);
  end
end

for i = 1:length(kern.comp)
  kern.nParams = kern.nParams + kern.comp{i}.nParams;
  kern.comp{i}.index = [];
end
kern.paramGroups = speye(kern.nParams);

if strcmp(kern.type, 'cmpnd')
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
end
