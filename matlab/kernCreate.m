function kern = kernCreate(X, kernelType)

% KERNCREATE Initialise a kernel structure.

% KERN

dim = size(X, 2);
if dim == 1 & size(X, 1) == 1;
  dim = X;
end
if iscell(kernelType)
  kern.inputDimension = dim;
  switch kernelType{1}
   case 'tensor'
    % tensor kernel type
    start = 2;
    kern.type = 'tensor';
   case 'cmpnd'
    % compound kernel type
    start = 2;
    kern.type = 'cmpnd';
   otherwise
    % compound kernel type
    start = 1;
    kern.type = 'cmpnd';
  end
  for i = start:length(kernelType)
    kern.comp{i-start+1} = kernCreate(X, kernelType{i});
    kern.comp{i-start+1}.index = [];
  end
  kern = kernParamInit(kern);
elseif isstruct(kernelType)
  % If a structure is passed, use it as the kernel.
  kern = kernelType;
else
  kern.type = kernelType;
  kern.inputDimension = dim;
  kern = kernParamInit(kern);
end
kern.Kstore = [];
kern.diagK = [];

