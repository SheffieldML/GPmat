function kern = kernCreate(X, kernelType)

% KERNCREATE Initialise a kernel structure.
% FORMAT
% DESC creates a kernel matrix structure given an design matrix of
% DESC input points and a kernel type.
% ARG X : Input data values (from which kernel will later be computed).
% ARG type : Type of kernel to be created, some standard types are
% 'lin', 'rbf', 'white', 'bias' and 'rbfard'. If a cell of the form
% {'cmpnd', 'rbf', 'lin', 'white'} is used a compound kernel
% based on the sum of the individual kernels will be created. The
% 'cmpnd' element at the start of the sequence is
% optional. Furthermore, {'tensor', 'rbf', 'lin'} can be used to
% give a tensor product kernel, whose elements are the formed from
% the products of the two indvidual kernel's elements and finally
% {'multi', 'rbf', ...} can be used to create a block structured
% kernel for use with multiple outputs.
% RETURN kern : The kernel structure.
%
% FORMAT
% DESC creates a kernel matrix structure given the dimensions of
% the design matrix and the kernel type.
% ARG dim : dimensions of the design matrix (in the form number of
% data points by number of features).
% ARG type : Type of kernel to be created, as described above.
% RETURN kern : The kernel structure.
%
% SEEALSO : kernParamInit

% KERN

dim = size(X, 2);
if dim == 1 & size(X, 1) == 1;
  dim = X;
end
if iscell(kernelType)
  kern.inputDimension = dim;
  switch kernelType{1}
   case 'multi'
    % multi output block based kernel.
    start = 2;
    kern.type = 'multi';
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

