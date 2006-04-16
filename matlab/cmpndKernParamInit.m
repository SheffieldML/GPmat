function kern = cmpndKernParamInit(kern)

% CMPNDKERNPARAMINIT Compound kernel parameter initialisation.
% FORMAT
% DESC initialises parameters on a compound kernel.
% ARG kern : the kernel structure to initialise.
% RETURN kern : the kernel structure returned with parameters
% initialised.
%
% SEEALSO : kernParamInit, tensorKernParamInit, cmpndKernCompute 

% KERN

kern.nParams = 0;
kern.transforms = [];
if ~isfield(kern, 'comp')
  kern.comp=cell(0);
end
for i = 1:length(kern.comp)
  kern.comp{i} = kernParamInit(kern.comp{i});
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

