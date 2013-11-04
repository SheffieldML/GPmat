function kern = pathKernCreate(X, symKernelType)

% PATHKERNCREATE Initialise a path kernel structure.
% FORMAT
% DESC creates a kernel matrix structure given an design matrix of
% DESC input points and a kernel type.
% ARG X : Input data values (from which kernel will later be computed).
% ARG symKernelType : type of symbol kernel, takes the same
% arguments as kernCreate
% 'lin', 'rbf', 'white', 'bias' and 'rbfard'. If a cell of the form
% {'cmpnd', 'rbf', 'lin', 'white'} is used a compound kernel
% based on the sum of the individual kernels will be created. The
% 'cmpnd' element at the start of the sequence is
% optional. Furthermore, {'tensor', 'rbf', 'lin'} can be used to
% give a tensor product kernel, whose elements are the formed from
% the products of the two indvidual kernel's elements and
% {'multi', 'rbf', ...} can be used to create a block structured
% kernel for use with multiple outputs.
% Finally the form {'parametric', struct('opt1', {val1}), 'rbf'} can
% be used to pass options to other kernels.
% RETURN kern : The kernel structure.
%
% SEEALSO : pathKernParamInit
%
% COPYRIGHT: Andrea Baisero, Carl Henrik Ek, 2013
%

% SHEFFIELDML


if(nargin<2)
    symKernelType='rbf';
end

kern.type='path';
kern.gkern=kernCreate(X,symKernelType);

kern.inputDimension=X;
kern=kernParamInit(kern);
