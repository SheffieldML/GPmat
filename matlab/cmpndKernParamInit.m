function kern = cmpndKernParamInit(kern)

% CMPNDKERNPARAMINIT Compound kernel parameter initialisation.

% IVM

kern.nParams = 0;
kern.transforms = [];
for i = 1:length(kern.comp)
  kern.comp{i} = kernParamInit(kern.comp{i});
  kern.nParams = kern.nParams + kern.comp{i}.nParams;
end
kern.paramGroups = speye(kern.nParams);

