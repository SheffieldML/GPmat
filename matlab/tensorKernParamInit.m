function kern = tensorKernParamInit(kern)

% TENSORKERNPARAMINIT Tensor kernel parameter initialisation.
% FORMAT
% DESC initialises parameters on a tensor product kernel.
% ARG kern : the kernel structure to initialise.
% RETURN kern : the kernel structure returned with parameters
% initialised.
%
% SEEALSO : kernParamInit, tensorKernCompute, cmpndKernParamInit 

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

% Warn if component kernels have white variance.
whiteVariance = 0;
for i = 1:length(kern.comp)
  if strcmp(kern.comp{i}.type, 'white')
    whiteVariance = whiteVariance + kern.comp{i}.variance;
  else
    if(isfield(kern.comp{i}, 'whiteVariance'))
      whiteVariance = whiteVariance + ...
          kern.comp{i}.whiteVariance;
    end
  end
end
if whiteVariance > 0
  warning(['Components of tensor kernel have non-zero white variance. ' ...
            'This can cause problems some implementations.']);
end
