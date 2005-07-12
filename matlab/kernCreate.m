function kern = kernCreate(X, kernelType)

% KERNCREATE Initialise a kernel structure.

% KERN

if iscell(kernelType)
  % compound kernel type
  kern.type = 'cmpnd';
  for i = 1:length(kernelType)
    kern.comp{i}.type = kernelType{i};
    kern.comp{i}.inputDimension = size(X, 2);
    kern.comp{i}.index = [];
  end
  kern = kernParamInit(kern);
elseif isstruct(kernelType)
  % If a structure is passed, use it as the kernel.
  kern = kernelType;
else
  kern.type = kernelType;
  kern.inputDimension = size(X, 2);
  kern = kernParamInit(kern);
end
kern.Kstore = [];
kern.diagK = [];

