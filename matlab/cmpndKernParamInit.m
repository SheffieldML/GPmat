function kern = cmpndKernParamInit(kern)

% CMPNDKERNPARAMINIT Compound kernel parameter initialisation.

% KERN

kern.nParams = 0;
kern.transforms = [];
for i = 1:length(kern.comp)
  kern.comp{i} = kernParamInit(kern.comp{i});
  kern.nParams = kern.nParams + kern.comp{i}.nParams;
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

