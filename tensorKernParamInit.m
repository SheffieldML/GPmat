function kern = tensorKernParamInit(kern)

% TENSORKERNPARAMINIT TENSOR kernel parameter initialisation.
% The tensor product (TENSOR) kernel is a container kernel for
% allowing tensor products kernels of separate component kernels.
% several different kernels to be added together. It is created by
% using the kernCreate command with the kernel type given as a
% cell and the first entry given as 'tensor'. For example, to
% create a tensor kernel compound kernel that is composed of
% an RBF kernel and a LIN kernel you call,
%
% kern = kernCreate(X, {'tensor', 'rbf', 'lin'});
%
% Each individual kernel is then stored within the returned kernel
% structure. The kernels are stored in order in a field called
% 'comp'. So display obtain the 'rbf' kernel you write:
%
% kernDisplay(kern.comp{1})
%
% SEEALSO : cmpndKernParamInit, multiKernParamInit
%
% FORMAT
% DESC initialises the tensor product
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006

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

% Warn if component kernels have white variance and check if
% resulting kernel is stationary.
kern.isStationary = true;
whiteVariance = 0;
for i = 1:length(kern.comp)
  if ~kern.comp{i}.isStationary
    kern.isStationary = false;
  end
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
